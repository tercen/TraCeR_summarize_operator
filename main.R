library(tercen)
library(dplyr)
library(readr)
library(stringr)

serialize.to.string = function(object){
  con = rawConnection(raw(0), "r+")
  saveRDS(object, con)
  str64 = base64enc::base64encode(rawConnectionValue(con))
  close(con)
  return(str64)
}
deserialize.from.string = function(str64){
  con = rawConnection(base64enc::base64decode(str64), "r+")
  object = readRDS(con)
  close(con)
  return(object)
}

find.schema.by.factor.name = function(ctx, factor.name){
  visit.relation = function(visitor, relation){
    if (inherits(relation,"SimpleRelation")){
      visitor(relation)
    } else if (inherits(relation,"CompositeRelation")){
      visit.relation(visitor, relation$mainRelation)
      lapply(relation$joinOperators, function(jop){
        visit.relation(visitor, jop$rightRelation)
      })
    } else if (inherits(relation,"WhereRelation") 
               || inherits(relation,"RenameRelation")){
      visit.relation(visitor, relation$relation)
    } else if (inherits(relation,"UnionRelation")){
      lapply(relation$relations, function(rel){
        visit.relation(visitor, rel)
      })
    } 
    invisible()
  }
  
  myenv = new.env()
  add.in.env = function(object){
    myenv[[toString(length(myenv)+1)]] = object$id
  }
  
  visit.relation(add.in.env, ctx$query$relation)
  
  schemas = lapply(as.list(myenv), function(id){
    ctx$client$tableSchemaService$get(id)
  })
  
  Find(function(schema){
    !is.null(Find(function(column) column$name == factor.name, schema$columns))
  }, schemas);
}


ctx <- tercenCtx()


schema <- find.schema.by.factor.name(ctx, names(ctx$cselect())[[1]])

table <- ctx$client$tableSchemaService$select(schema$id, Map(function(x) x$name, schema$columns), 0, schema$nRows)

table <- as_tibble(table)

hidden_colnames <- colnames(table)[str_starts(colnames(table), "\\.")]


if (!((".forward_read_fastq_data" %in% hidden_colnames) &
      (".reverse_read_fastq_data" %in% hidden_colnames))) {
  
  stop("Input is not samples containing paired-end fastq data.")
  
}

for (i in 1:nrow(table)) {
  
  sample_name <- select(table, ends_with(".sample"))[[1]][[i]]
  
  filename_r1 <- paste0(sample_name, "1.fastq.gz")
  filename_r2 <- paste0(sample_name, "2.fastq.gz")
  
  writeBin(deserialize.from.string(table[".forward_read_fastq_data"][[1]][[i]]), filename_r1)
  writeBin(deserialize.from.string(table[".reverse_read_fastq_data"][[1]][[i]]), filename_r2)
  
  cmd = '/tracer/tracer'
  args = paste('assemble',
               '--ncores', parallel::detectCores(),
               '--config_file /tercen_tracer.conf',
               '-s Hsap',
               filename_r1, filename_r2,
               paste0(sample_name, "_", i), "this_run",
               sep = ' ')
  
  system2(cmd, args)
  
}

args <- paste('summarise',
              '--ncores', parallel::detectCores(),
              '--config_file /tercen_tracer.conf',
              '-s Hsap',
              'this_run',
              sep = ' ')

system2(cmd, args)

# run the collect script

cmd <- "python /collect_TRA_TRB_in_fasta.py this_run/*/filtered_TCR_seqs/*.fa > this_run.tsv"

system(cmd)

collected_summary <- read_tsv("this_run.tsv")

cols <- sapply(collected_summary, is.logical)
collected_summary[,cols] <- lapply(collected_summary[,cols], as.numeric)


(collected_summary %>%
    mutate(.ci = 0) %>%
    ctx$addNamespace() %>%
    ctx$save())

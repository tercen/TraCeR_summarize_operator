library(tercen)
library(dplyr)
library(readr)

ctx <- tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required") 

documentIds <- ctx$cselect("documentId")

file_names <- sapply(documentIds[[1]],
                     function(x) (ctx$client$fileService$get(x))$name) %>%
  sort()

if((length(file_names) %% 2) != 0) stop("Non-even number of files supplied. Are you sure you've supplied paired-end files?")


for (first_in_pair_index in seq(1, length(file_names), by = 2)) {
  
  docIds = file_names[first_in_pair_index:(first_in_pair_index+1)]
  
  dodId_r1 <- names(docIds)[[1]]
  doc_r1 <- ctx$client$fileService$get(dodId_r1)
  filename_r1 <- docIds[[1]]
  writeBin(ctx$client$fileService$download(dodId_r1), filename_r1)
  on.exit(unlink(filename_r1))
  
  dodId_r2 <- names(docIds)[[2]]
  doc_r2 <- ctx$client$fileService$get(dodId_r2)
  filename_r2 <- docIds[[2]]
  writeBin(ctx$client$fileService$download(dodId_r2), filename_r2)
  on.exit(unlink(filename_r2))
  
  seqlibrary_name <- substr(filename_r1, 1,which.min(strsplit(filename_r1, "")[[1]] == strsplit(filename_r2, "")[[1]]))
  
  cmd = '/tracer/tracer'
  args = paste('assemble',
               '--ncores', parallel::detectCores(),
               '--config_file /tercen_tracer.conf',
               '-s Hsap',
               filename_r1, filename_r2,
               seqlibrary_name, "this_run",
               sep = ' ')
  
  system2(cmd, args)

  stop("Hammertime")
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

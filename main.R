library(tercen)
library(dplyr)
library(readr)
library(stringr)


ctx <- tercenCtx()



#if (length(ctx$cselect()) == 2 ) {
#  sample_names_vector = ctx$cselect()[[2]]
#} else {
#  sample_names_vector = ctx$cselect()[[1]]
#}

# Define input path
# LATER THIS WILL COME FROM THE TERCEN TABLE!
input_path <- "/var/lib/tercen/external/read/admin/tracer_summarise_dev"

print(input_path)

# Define output path
# LATER THIS WILL COME FROM THE TERCEN TABLE!
output_path <- "/var/lib/tercen/external/read/admin/tracer_summarise_dev"

# run the TraCeR summarise command

cmd = '/tracer/tracer'

args <- paste('summarise',
              '--ncores', parallel::detectCores(),
              '--config_file /tercen_tracer.conf',
              '-s Hsap',
              input_path,
              sep = ' ')

system2(cmd, args)

# run the collect script

cmd <- paste0("python /collect_TRA_TRB_in_fasta.py ", input_path, "/*/filtered_TCR_seqs/*.fa > ", output_path, "/summarise_output.tsv")

system(cmd)

# move the tracer_summarise output to the output folder
cmd <- paste0("mv ", input_path, "/filtered_TCRAB_summary ", output_path, "/")

system(cmd)

collected_summary <- read_tsv(paste0(input_path, "/summarise_output.tsv"))

recombinants <- read_tsv(paste0(input_path, "/filtered_TCRAB_summary/recombinants.txt")

collected_summary <- left_join(collected_summary,
                               recombinants,
                               by = c(sample = "cell_name"))

cols <- sapply(collected_summary, is.logical)
collected_summary[,cols] <- lapply(collected_summary[,cols], as.numeric)

#save_output <- as.character(ctx$op.value('save_output_to_folder'))
#
#if (save_output == "yes") {
#  
#  output_folder_prefix <- as.character(ctx$op.value('output_folder_prefix'))
#  
#  # create trim galore zipped output
#  system("zip -r tracer_output.zip this_run")
#  
#  # save zipped file to project folder
#  filename <- "tracer_output.zip"
#  bytes = readBin(file(filename, 'rb'),
#                  raw(),
#                  n = file.info(filename)$size)
#  
#  fileDoc = FileDocument$new()
#  fileDoc$name = paste0(output_folder_prefix, "_", filename)
#  fileDoc$projectId = ctx$cschema$projectId
#  fileDoc$size = length(bytes)
#  
#  fileDoc = ctx$client$fileService$upload(fileDoc, bytes)
#  
#}


(collected_summary %>%
    mutate(.ci = 0) %>%
    ctx$addNamespace() %>%
    ctx$save())

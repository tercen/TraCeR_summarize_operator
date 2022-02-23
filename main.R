library(tercen)
library(dplyr)
library(readr)
library(stringr)


ctx <- tercenCtx()


# Define output path
output_path <- "/var/lib/tercen/share/write/tracer_output"

# Check if a "tracer_output" folder exists and is not empty
if( dir.exists(output_path) == FALSE) {

  stop("ERROR: tracer_output folder does not exist in project write folder.")

}

if (length(dir(output_path)) == 0) {
  stop("ERROR: tracer_output folder is empty.")
}

# run the TraCeR summarise command

cmd = '/tracer/tracer'

args <- paste('summarise',
              '--ncores', parallel::detectCores(),
              '--config_file /tercen_tracer.conf',
              '-s Hsap',
              output_path,
              sep = ' ')

system2(cmd, args)

# run the collect script

cmd <- paste0("python /collect_TRA_TRB_in_fasta.py ", output_path, "/*/filtered_TCR_seqs/*.fa > ", output_path, "/summarise_output.tsv")

system(cmd)

collected_summary <- read_tsv(paste0(output_path, "/summarise_output.tsv"))

recombinants <- read_tsv(paste0(output_path, "/filtered_TCRAB_summary/recombinants.txt"))

collected_summary <- left_join(collected_summary,
                               recombinants,
                               by = c(sample = "cell_name"))

cols <- sapply(collected_summary, is.logical)
collected_summary[,cols] <- lapply(collected_summary[,cols], as.numeric)


(collected_summary %>%
    mutate(.ci = 0) %>%
    ctx$addNamespace() %>%
    ctx$save())

library(tercen)
library(dplyr)

ctx = tercenCtx()

if (!any(ctx$cnames == "documentId")) stop("Column factor documentId is required")

df <- ctx$cselect()

command_to_run <- "/tracer/tracer test -c tercen_tracer.conf"

system(command_to_run)

tibble(.ci = c(0),
       status = c(1)) %>%
ctx$addNamespace() %>%
ctx$save()

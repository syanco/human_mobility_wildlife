if(interactive()) {
  library(here)
  
  .wd <- '~/projects/covid-19_movement'
  .test <- TRUE
  # rd <- here::here
  
  .datPF <- file.path(.wd,'analysis/src/workflow/')
  .outPF <- file.path(.wd,"analysis/src/hpc/")
  
} else {
  library(docopt)
  # library(rprojroot)
  
  .wd <- '/gpfs/loomis/pi/jetz/sy522/covid-19_movement'
  # .script <-  thisfile()
  # rd <- is_rstudio_project$make_fix_file(.script)
  
  .datPF <- file.path(.wd,'analysis/src/workflow/')
  .outPF <- file.path(.wd,"analysis/src/hpc/")
}

source(file.path(.wd,'analysis/src/startup.r'))

suppressWarnings(
  suppressPackageStartupMessages({
    library(data.table)
  }))



files <- list.files(paste0(.wd,'/out/ssf-background-pts/individual-files'),full.names = TRUE)

n_total <- length(files)


n <- 150
n_events <- round(n_total/n)
# n_events <- ifelse((n_total/n)==n_events, n_events, n_events)

start_ix <- seq(from = 1, to = n_total, by = n_events)
end_ix <- c(seq(from = n_events, to = n_total, by = n_events),n_total)


joblist <- data.frame("string" = rep(paste0(" module load miniconda; conda activate move; Rscript ",.datPF,"annotate-background-sg-ghm.r "), times = n),
                      "arg1" = start_ix,
                      "arg2" = end_ix,
                      "arg3" = seq(from = 1, to = n, by = 1))

write.table(joblist,paste0(.outPF,"annotation-joblist.txt"),col.names = FALSE,row.names = FALSE,quote = FALSE)
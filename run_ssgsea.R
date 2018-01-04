sessionInfo()
options(echo=FALSE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
exp_file <- args[1]
gmt_file <- args[2]
weight <- as.numeric(args[3])
out_file <- args[4]
ccba.path <- args[5]
source(ccba.path)
CCBA_ssGSEA_project_dataset.v1(input.ds=exp_file, output.ds=out_file, gene.set.databases=gmt_file, weight=weight)
  
library(oligo)

# Get the input directory and output file from the command line
args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_file <- args[2]

# Read the CEL files from the specified input directory without normalization
celfiles <- list.files(input_dir, pattern = "*.CEL.gz", full.names=TRUE)
raw_data <- read.celfiles(celfiles)

rma_normalized <- rma(raw_data)

# Write the RMA normalized expression values to the output file
write.exprs(rma_normalized, file=output_file)

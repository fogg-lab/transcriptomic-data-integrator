library(affy)

# Get the input directory and output file from the command line
args <- commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_file <- args[2]

# Read the CEL files from the specified input directory
data <- ReadAffy(celfile.path=input_dir)

# Perform RMA normalization
rma_normalized <- rma(data)

# Write the RMA normalized expression values to the output file
write.exprs(rma_normalized, file=output_file)

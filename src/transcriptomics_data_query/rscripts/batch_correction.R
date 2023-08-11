library(readr)
library(sva)

# Read in the arguments
args <- commandArgs(trailingOnly = TRUE)
data_type <- tolower(args[1])
expression_matrix_file <- args[2]
clinical_data_file <- args[3]
output_file <- args[4]

# Verify that the data type is valid (either "rnaseq" or "microarray")
if(!(data_type %in% c("rnaseq", "microarray"))) {
  stop("Invalid data type. Must be either 'rnaseq' or 'microarray'.")
}

# Read in the sample data and expression
expression <- read_tsv(expression_matrix_file)
clinical_data <- read_tsv(clinical_data_file)

# Prep batch correction
rma_expr <- as.matrix(expression[-1])
rownames(rma_expr) <- expression$symbol
model_m <- model.matrix(~ condition, data = clinical_data)
batch <- clinical_data$batch

# Get symbol column for adding to output later
symbols <- data.frame(symbol = expression$symbol)

# Batch correction
if (data_type == "rnaseq") {
  bc_rma_expr <- ComBat_seq(rma_expr, batch = batch, group = NULL, covar_mod=model_m)
} else { # microarray
  bc_rma_expr <- ComBat(rma_expr, batch = batch, mod = model_m, ref.batch = 1)
}

# Add symbol column and write the dataframe to a file
bc_rma_expr <- data.frame(bc_rma_expr)
result <- dplyr::bind_cols(symbols, bc_rma_expr)
write_tsv(result, file.path(output_dir, "batch_corrected_expression.tsv"))

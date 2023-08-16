library(readr)
library(sva)

# Read in the arguments
args <- commandArgs(trailingOnly = TRUE)
data_type <- tolower(args[1])
expression_matrix_file <- args[2]
clinical_data_file <- args[3]
condition_column <- args[4]
output_file <- args[5]

# Verify that the data type is valid (either "rnaseq" or "microarray")
if(!(data_type %in% c("rnaseq", "microarray"))) {
  stop("Invalid data type. Must be either 'rnaseq' or 'microarray'.")
}

# Read in the sample data and expression
expression_df <- read_tsv(expression_matrix_file)
clinical_data <- read_tsv(clinical_data_file)

# Check if the specified condition column is in the clinical data table
if (!(condition_column %in% names(clinical_data))) {
  stop(paste("Specified condition column '", condition_column, "' not found in clinical data table."))
}

# Prep batch correction
expression_matrix <- as.matrix(expression_df[-1])
rownames(expression_matrix) <- expression_df[[1]]
model_m <- model.matrix(as.formula(paste("~", condition_column)), data = clinical_data)
batch <- clinical_data$batch

# Get the first column for adding to output later
symbols <- data.frame(symbol = expression_df[[1]])

# Batch correction
if (data_type == "rnaseq") {
  bc_expression_matrix <- ComBat_seq(expression_matrix, batch = batch, group = NULL, covar_mod=model_m)
} else { # microarray
  bc_expression_matrix <- ComBat(expression_matrix, batch = batch, mod = model_m, ref.batch = 1)
}

# Add the first column and write the dataframe to a file
bc_expression_matrix <- data.frame(bc_expression_matrix)
result <- dplyr::bind_cols(symbols, bc_expression_matrix)
write_tsv(result, output_file)

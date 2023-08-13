library(readr)

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

normalization_method <- tolower(args[1])
expression_matrix_file <- args[2]
clinical_data_file <- args[3]
design_column <- args[4]
output_file <- args[5]

# Verify that the normalization method is valid (either "mrn" or "tmm")
if(!(normalization_method %in% c("mrn", "tmm"))) {
    stop("Invalid normalization method. Must be either 'mrn' or 'tmm'.")
}

expression <- read_tsv(expression_matrix_file, show_col_types = FALSE)
clinical_data <- read_tsv(clinical_data_file, show_col_types = FALSE)

# Verify that the design column is present in the clinical data
if (!(design_column %in% names(clinical_data))) {
  stop(paste("Design column", design_column, "not found in clinical data."))
}

# Select only numeric columns
numeric_expression <- expression[, sapply(expression, is.numeric)]

expr_mat <- as.matrix(numeric_expression)

# initialize empty matrix as a placeholder for the normalized expression
normalized_expression <- matrix(, nrow = 1, ncol = 1)

if(normalization_method == "mrn") {
  # Verify that the design column is present in the clinical data
  if (!(design_column %in% names(clinical_data))) {
    stop(paste("Design column", design_column, "not found in clinical data."))
  }

  library(DESeq2)
  design_formula <- as.formula(paste("~", design_column))
  dds <- DESeqDataSetFromMatrix(countData = expr_mat, colData=clinical_data, design = design_formula)
  dds <- estimateSizeFactors(dds)
  normalized_expression <- expression(dds, normalized=TRUE)
} else if (normalization_method == "tmm") {
  library(edgeR)
  dge_list <- DGEList(expression=expr_mat, samples=clinical_data)
  dge_list <- calcNormFactors(dge_list, method="TMM")
  normalized_expression <- cpm(dge_list)
}

normalized_expression <- as.data.frame(normalized_expression)
write.table(normalized_expression, file=output_file, sep="\t", quote=F, row.names=FALSE)

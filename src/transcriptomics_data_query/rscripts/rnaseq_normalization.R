library(readr)
library(edgeR)

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

expression_matrix_file <- args[1]
clinical_data_file <- args[2]
output_file <- args[3]

expression <- read_tsv(expression_matrix_file, show_col_types = FALSE)
clinical_data <- read_tsv(clinical_data_file, show_col_types = FALSE)

expr_mat <- as.matrix(expression[, sapply(expression, is.numeric)])

dge_list <- DGEList(counts=expr_mat, samples=clinical_data)
dge_list <- calcNormFactors(dge_list, method="TMM")
normalized_expression <- as.data.frame(cpm(dge_list))

write.table(normalized_expression, file=output_file, sep="\t", quote=F, row.names=FALSE)

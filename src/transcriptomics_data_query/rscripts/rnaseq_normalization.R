library(readr)

# command line arguments
args <- commandArgs(trailingOnly = TRUE)

normalization_method <- tolower(args[1])
expression_matrix_file <- args[2]
clinical_data_file <- args[3]
output_file <- args[4]

# Verify that the normalization method is valid (either "mrn" or "tmm")
if(!(normalization_method %in% c("mrn", "tmm"))) {
    stop("Invalid normalization method. Must be either 'mrn' or 'tmm'.")
}

expression <- read_tsv(expression_matrix_file, show_col_types = FALSE)
clinical_data <- read_tsv(clinical_data_file, show_col_types = FALSE)
genes <- expression[,1, drop=FALSE]

# filter out unnecessary columns
cols = colnames(expression)
extraneous_cols = c("hugo_symbol", "entrez_gene_id", "symbol", "gene")
for (x in 1:length(cols)) {
    if(tolower(cols[x]) %in% extraneous_cols) {
        expression <- expression[, -which(names(expression) == cols[x])]
    }
}

expr_mat <- as.matrix(expression)

# initialize empty matrix as a placeholder for the normalized expression
normalized_expression <- matrix(, nrow = 1, ncol = 1)

if(normalization_method == "mrn") {
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = expr_mat, colData=clinical_data, design = ~condition)
    dds <- estimateSizeFactors(dds)
    normalized_expression <- expression(dds, normalized=TRUE)
} else if (normalization_method == "tmm") {
    library(edgeR)
    dge_list <- DGEList(expression=expr_mat, samples=clinical_data)
    dge_list <- calcNormFactors(dge_list, method="TMM")
    normalized_expression <- cpm(dge_list)
}

rownames(normalized_expression) <- NULL
normalized_expression <- as.data.frame(normalized_expression)
final_expression <- do.call(cbind, list(genes, normalized_expression))
write.table(final_expression, file=output_file, sep="\t", quote=F, row.names=FALSE)

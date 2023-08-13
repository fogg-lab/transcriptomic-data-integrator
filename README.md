# transcriptomics-data-query-and-retrieval
Search, retrieve, and prepare transcriptomics data from GEO and GDC (development in progress).

## Prerequisites
- Python 3.7 or higher.
- To use the functions for normalization or batch correction from the `preprocess` module, install [R](https://www.r-project.org/) and packages:
  - [affy](https://bioconductor.org/packages/release/bioc/html/affy.html)
  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [readr](https://cran.r-project.org/web/packages/readr/index.html)
  - [sva](https://bioconductor.org/packages/release/bioc/html/sva.html)
  - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

These packages can be installed in an R environment using the following commands:
```R
install.packages('BiocManager', repos='https://cloud.r-project.org/')
BiocManager::install(c('affy', 'edgeR', 'readr', 'sva', 'DESeq2'))
```

## Setup

Run the below commands at the command line. Replace dummy email with your email which will be submitted in your GEO queries to the NCBI API.
```zsh
git clone https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval.git
cd transcriptomics-data-query-and-retrieval
pip install -e .
configure-ncbi-email YOUR_EMAIL@EXAMPLE.COM
```

## Examples

In progress: Colab notebook with examples.

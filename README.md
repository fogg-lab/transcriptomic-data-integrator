# transcriptomics-data-query
A Python package to search, retrieve, and prepare gene expression data from Gene Expression Omnibus and Genomic Data Commons.

The `gdc` module (to query/retrieve data from GDC) is not implemented yet. It is currently in development.

## Prerequisites
- Python 3.7 or higher.
- To use the functions for normalization or batch correction from the `preprocess` module, install [R](https://www.r-project.org/) and packages:
  - [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
  - [oligo](https://www.bioconductor.org/packages/release/bioc/html/oligo.html)
  - [readr](https://cran.r-project.org/web/packages/readr/index.html)
  - [sva](https://bioconductor.org/packages/release/bioc/html/sva.html)

For RMA normalization, you will need to install platform design info packages, such as:
  - [pd.clariom.d.human](https://bioconductor.org/packages/release/data/annotation/html/pd.clariom.d.human.html)
  - [pd.hg.u133.plus.2](https://bioconductor.org/packages/release/data/annotation/html/pd.hg.u133.plus.2.html)
  - other packages for different platforms you might encounter

These packages can be installed in an R environment by running the script install_r_packages.R. This install script was written for R 4.3.

## Setup: Install and configure the `transcriptomics_data_query` package

Run the below commands at the command line. Replace dummy email with your email which will be submitted in your GEO queries to the NCBI API.
```zsh
git clone https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval.git
cd transcriptomics-data-query-and-retrieval
pip install -e .
configure-ncbi-email YOUR_EMAIL@EXAMPLE.COM
```

## Usage

Refer to the [documentation](https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/blob/main/DOCUMENTATION.md) and [Colab notebooks](https://github.com/fogg-lab/transcriptomics-data-query-and-retrieval/tree/main/notebooks).

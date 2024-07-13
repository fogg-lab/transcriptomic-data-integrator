# transcriptomic-data-integrator
A Python package to retrieve and prepare gene expression data from Gene Expression Omnibus and Genomic Data Commons.

**Archived**: This was my earlier attempt to automate common data retrieval and preprocessing tasks for gene expression data. Interested readers might want to check out these other resources instead:
- https://github.com/alexvpickering/crossmeta
- https://maayanlab.cloud/archs4/
- https://jhubiostatistics.shinyapps.io/recount/

The `gdc` module (to retrieve and prepare data from GDC) is not implemented.

## Prerequisites
- Python 3.9 or higher.
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

## Setup: Install and configure the `transcriptomic_data_integrator` package

Run the below commands at the command line. Replace dummy email with your email which will be submitted in your GEO queries to the NCBI API.
```zsh
git clone https://github.com/fogg-lab/transcriptomic-data-integrator.git
cd transcriptomic-data-integrator
pip install -e .
configure-ncbi-email YOUR_EMAIL@EXAMPLE.COM
```

## Usage

Refer to the [documentation](https://github.com/fogg-lab/transcriptomic-data-integrator/blob/main/DOCUMENTATION.md) and [Colab notebooks](https://github.com/fogg-lab/transcriptomic-data-integrator/tree/main/notebooks).

## Known limitations

- The function `tdi.geo.map_probes_to_genes` is not guaranteed to work on all microarray platform technologies. This is due to differences in how the probe set annotation table is organized between different platforms.
- Other GEO query functions, such as `tdi.geo.get_geo_clinical_characteristics`, fail when the data for the study on GEO is not organized according to how this package expects. This happens more times than not.

If you encounter any problems using the package, please [submit an issue](https://github.com/fogg-lab/transcriptomic-data-integrator/issues/new) to report it.

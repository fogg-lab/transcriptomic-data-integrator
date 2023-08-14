install.packages(c('BiocManager', 'readr'), repos='https://cloud.r-project.org/', Ncpus = 2)
BiocManager::install(version = "3.17")
BiocManager::install(c('affy', 'edgeR', 'sva', 'DESeq2', 'oligo', 'pd.clariom.d.human'), Ncpus = 2)

install.packages(c('BiocManager', 'readr'), repos='https://cloud.r-project.org/', Ncpus = 2)
BiocManager::install(version = "3.17")
BiocManager::install(c('edgeR', 'sva', 'oligo', 'pd.clariom.d.human', 'pd.hg.u133.plus.2'), Ncpus = 2)

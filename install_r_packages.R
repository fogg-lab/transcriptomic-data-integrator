install.packages(c('BiocManager', 'readr'), repos='https://cloud.r-project.org/', Ncpus = 2)

# Install preprocessCore without threading to avoid a multithreading bug using the rma function
BiocManager::install('preprocessCore', configure.args='--disable-threading', force = TRUE)

BiocManager::install(c('affy', 'edgeR', 'sva', 'DESeq2', 'oligo', 'pd.clariom.d.human', Ncpus = 2))

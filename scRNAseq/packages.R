# Use this script to install packages via CRAN, for example: 

# NOTE: Code below shows how you install R packages from CRAN and Bioconductor. For CRAN packages, you can use 
# the standard install.packages() function; for Bioconductor packages, however, you 
# must first install BiocManager and then use that for installs. 

#update.packages(repos='http://cran.us.r-project.org', ask=FALSE)

#install.packages('remotes', dependencies = TRUE, repos = repos='http://cran.rstudio.com/')
#library('remotes')
#remotes::install_version(package = 'Seurat', version = package_version('3.2.1'))
#library('devtools')
#devtools::install_github("immunogenomics/harmony")
#install.packages("patchwork", repos='http://cran.rstudio.com/')
#install.packages("tidyr", repos='http://cran.rstudio.com/')
#install.packages("ggplot2", repos='http://cran.rstudio.com/')
BiocManager::install("AnnotationHub", ask = FALSE, update=FALSE)
BiocManager::install("biomaRt", update = FALSE, ask = FALSE)
#install.packages("ggplot2", repos='http://cran.rstudio.com/')
#install.packages('sctransform', repos='http://cran.rstudio.com/')

install.packages(c("broom", "corrplot", "car", "cowplot", "MASS", "fitdistrplus", "future", "future.apply", "ggrepel", "ggridges", "httr", "ica", "igraph", "irlba", "jsonlite", "leiden", "lmtest", "matrixStats", "miniUI", "patchwork", "pbapply", "plotly", "png", "RANN", "Rcpp", "RcppAnnoy", "reticulate", "ROCR", "rsvd", "Rtsne", "sctransform", "shiny", "spatstat", "uwot", "RcppEigen", "RcppProgress"),  repos = 'https://cloud.r-project.org')

install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos=NULL,type="source")
install.packages('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.2.1.tar.gz', repos = NULL, type = "source")
install.packages('https://cran.r-project.org/src/contrib/Archive/rstatix/rstatix_0.6.0.tar.gz', repos = NULL, type = "source")
install.packages('https://cran.r-project.org/src/contrib/Archive/sctransform/sctransform_0.2.1.tar.gz', repos = NULL, type = "source")
#install.packages('https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.7.tar.gz', repos = NULL, type = "source")


FROM rocker/rstudio:4.0.2

USER root

WORKDIR /home/rstudio

RUN apt-get update && \
	apt-get -y --no-install-recommends install --fix-missing \
        apt-file \
        apt-utils \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
	    default-jdk \
        default-jre\
        gdb \
        git \
        hdf5-helpers \
        lbzip2 \
        libbz2-dev \
        libcairo2-dev \
        libcurl4-openssl-dev \
        libfftw3-dev \
        libgeos-dev \
        libgl1-mesa-dev \
        libglpk-dev \
        libglu1-mesa-dev \
        libgsl0-dev \
        libhdf4-alt-dev \
        libhdf5-dev \
        libjpeg-dev \
        libjq-dev \
        liblzma-dev \
        libmariadbd-dev \
        libnetcdf-dev \
        libpng-dev \
        libpq-dev \
        libproj-dev \
        libprotobuf-dev \
        libsqlite3-dev \
        libssh2-1-dev \
        libssl-dev \
        libudunits2-dev \
        libxml2-dev \
        libxt-dev \
        libz-dev \
        make \
        netcdf-bin \
        postgis \
        protobuf-compiler \
        python3-pip \
        sqlite3 \
        tk-dev \
        unixodbc-dev \
        unzip \
        vim \
        libpoppler-cpp-dev \
        && apt-get clean && rm -rf /var/lib/apt/lists/*

# FFTW, FIt-SNE
RUN wget http://www.fftw.org/fftw-3.3.8.tar.gz \

    && tar zxvf fftw-3.3.8.tar.gz \
    && rm fftw-3.3.8.tar.gz \
    && cd fftw-3.3.8 \
    && ./configure \
    && make \
    && make install \
    && git clone https://github.com/KlugerLab/FIt-SNE.git \
    && cd FIt-SNE/ \
    && g++ -std=c++11 -O3  src/sptree.cpp src/tsne.cpp src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm \
    && cp bin/fast_tsne /usr/local/bin/

# Install R packages
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/', ask = FALSE)"
RUN R -e "install.packages('devtools', dependencies=TRUE, repos='http://cran.rstudio.com/', ask = FALSE)"

RUN R -e 'install.packages(c("broom", "corrplot", "car", "cowplot", "MASS", "fitdistrplus", "future", "future.apply", "ggrepel", "ggridges", "hdf5r", "httr", "ica", "igraph", "irlba", "jsonlite", "leiden", "lmtest", "matrixStats", "miniUI", "patchwork", "pbapply", "plotly", "png", "RANN", "Rcpp", "RcppAnnoy", "reticulate", "ROCR", "rsvd", "Rtsne", "sctransform", "shiny", "spatstat", "tidyverse", "uwot", "RcppEigen", "RcppProgress", "RColorBrewer", "viridis", "spatstat", "rstatix"), dependencies=TRUE, repos="http://cran.rstudio.com/", ask = FALSE)' \

RUN R -e 'BiocManager::install(c( "minfi", "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "batchelor", "AnnotationHub", "biomaRt", "clusterProfiler", "enrichplot", "batchelor", "Biobase", "limma", "rhdf5filters", "XVector", "GenomeInfoDb", "sparseMatrixStats", "BiocNeighbors", "DelayedArray", "rhdf5", "DOSE", "beachmat", "GenomicRanges", "ResidualMatrix", "HDF5Array", "BiocSingular", "scuttle"),  update = FALSE, ask = FALSE)' \

RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz"", repos=NULL,type="source")' \
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.2.1.tar.gz", repos = NULL, type = "source")' \
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/rstatix/rstatix_0.6.0.tar.gz", repos = NULL, type = "source")' \
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/sctransform/sctransform_0.2.1.tar.gz", repos = NULL, type = "source")' \

RUN R -e 'devtools::install_github("cole-trapnell-lab/leidenbase")' \
RUN R -e 'devtools::install_github("cole-trapnell-lab/monocle3")' 

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
        r-cran-rjava \ 
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

RUN R -e "install.packages(c('broom', 'corrplot', 'car', 'cowplot', 'MASS', 'fitdistrplus', 'future', 'future.apply', 'ggrepel', 'ggridges', 'httr', 'ica', 'igraph', 'irlba', 'jsonlite', 'leiden', 'lmtest','knitr', 'kableExtra','matrixStats', 'miniUI', 'patchwork', 'pbapply', 'plotly', 'png', 'RANN', 'Rcpp', 'RcppAnnoy', 'reticulate', 'ROCR', 'rsvd', 'Rtsne', 'sctransform', 'tidyverse', 'RColorBrewer', 'viridis', 'spatstat', 'uwot', 'RcppProgress'), dependencies=TRUE, repos='http://cran.rstudio.com/', ask = FALSE)"

RUN R -e "BiocManager::install(c('AnnotationHub','DOSE','GenomeInfoDb','GenomicRanges','SingleCellExperiment','SummarizedExperiment','biomaRt','clusterProfiler','enrichplot', 'limma', 'edgeR'), update = FALSE)"

RUN R -e "devtools::install_github('colearendt/xlsx')"

RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/rstatix/rstatix_0.6.0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/crayon/crayon_1.3.4.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/hdf5r/hdf5r_1.3.2.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/R6/R6_2.4.1.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/SeuratObject/SeuratObject_4.0.0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/stringi/stringi_1.5.3.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/withr/withr_2.3.0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.3-0.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/spatstat/spatstat_1.64-1.tar.gz', repos = NULL, type = 'source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.2.2.tar.gz', repos = NULL, type = 'source')"

#COPY ./packages.R /home/rstudio

#RUN Rscript /home/rstudio/packages.R

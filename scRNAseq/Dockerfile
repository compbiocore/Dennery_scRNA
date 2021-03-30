FROM rocker/rstudio:4.0.2
#FROM jupyter/r-notebook:54462805efcb

USER root

WORKDIR /home/rstudio

#install some base libs
RUN apt-get update && \
	apt-get -y --no-install-recommends install --fix-missing \
	gdb \
	vim \
	libxml2-dev \
	libcairo2-dev \
	libsqlite3-dev \
	libmariadbd-dev \
	libpq-dev \
	libssh2-1-dev \
	unixodbc-dev \
	libcurl4-openssl-dev \
	libssl-dev \
	libbz2-dev \
	liblzma-dev \
	python3-pip \
	libz-dev \
	libpng-dev \
        default-jre\
	default-jdk \
	libjpeg-dev \
	libglpk-dev \
	libxt-dev \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# Add BiocVersion
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/', ask = FALSE)" \
        && R -e "BiocManager::install(version='3.11', ask = FALSE)" 

COPY ./packages.R /home/rstudio
#COPY ./Seurat_3.2.1.tar.gz /home/rstudio
#COPY ./Analysis.Rmd /home/rstudio
#COPY ./MCA /home/MCA
#COPY ./O2D60out /home/O2D60out
#COPY ./O2D7out /home/O2D7out
#COPY ./AirD7out /home/AirD7out
#COPY ./AirD60out /home/AirD60out

RUN Rscript /home/rstudio/packages.R
#ADD https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.2.1.tar.gz
#R CMD INSTALL Seurat_3.2.1.tar.gz
#RUN R -e "install.packages("Seurat_3.2.1.tar.gz", repos = NULL, type = "source")




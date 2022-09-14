
## BIOCONDUCTOR
## Pull the development version of the Bioconductor Docker image
FROM bioconductor/bioconductor_docker:devel

RUN apt-get update \
	## Install the python package sceptre (for schoof2021 vignette)
	&& pip install git+https://github.com/bfurtwa/Sceptre.git@818a8914fe87788642f9b0dcdb49991ba8a4506a \
    ## Install specific version of NumPy to solve dependency issues, and install 
	## other python dependencies
    && pip install NumPy==1.22 IPython leidenalg \ 
	## Remove packages in '/var/cache/' and 'var/lib'
	## to remove side-effects of apt-get update
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/* \
	## Switch to libblas for better integration of python in reticulate. 
	&& ARCH=$(uname -m) \
	&& update-alternatives --set "libblas.so.3-${ARCH}-linux-gnu" "/usr/lib/${ARCH}-linux-gnu/blas/libblas.so.3" \
	&& update-alternatives --set "liblapack.so.3-${ARCH}-linux-gnu" "/usr/lib/${ARCH}-linux-gnu/lapack/liblapack.so.3"

## Copy the SCP.replication repo to give access to the vignette source code
COPY . /home/rstudio/SCP.replication/

## Change home directory
WORKDIR /home/rstudio/SCP.replication/vignettes/

## Install SCP.replication and dependencies, and copy source 
RUN R -e "BiocManager::install('UCLouvain-CBIO/SCP.replication', dependencies = TRUE)"

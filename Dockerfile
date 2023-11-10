
## BIOCONDUCTOR
## Pull the development version of the Bioconductor Docker image
FROM bioconductor/bioconductor_docker:devel

## Change home directory
WORKDIR /home/rstudio/

## Install QFeatures, scpdata and their dependencies (2 passes)
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
RUN R -e "BiocManager::install(c('rformassspectrometry/QFeatures', 'UCLouvain-CBIO/scpdata'), dependencies = TRUE)"
## Install latest scp version (containing the scplainer implementation) and dependencies (2 passes)
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', ref = 'asca_merge', dependencies = TRUE)"
RUN R -e "BiocManager::install('UCLouvain-CBIO/scp', ref = 'asca_merge', dependencies = TRUE)"
## Install SCP.replication and dependencies (2 passes)
RUN R -e "BiocManager::install('UCLouvain-CBIO/SCP.replication', dependencies = TRUE)"
RUN R -e "BiocManager::install('UCLouvain-CBIO/SCP.replication', dependencies = TRUE)"

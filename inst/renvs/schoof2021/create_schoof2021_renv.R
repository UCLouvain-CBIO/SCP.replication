####---- Initialize renv ----####

## In inst/renvs/schoof2021 directory
## It should be only run once
renv::init(bare = TRUE)

####---- Activate renv ----####

## Run this every time you want to use the renv or to modify it
renv::activate("~/PhD/SCP.replication/inst/renvs/schoof2021/")

## In order to include a package in an `renv`, it needs:
## 1. to be installed locally
## 2. to be loaded from a script using the `library` function

####---- Install packages ----####

renv::install("BiocManager")
BiocManager::install(version = "devel")
## Add Bioconductor repos. By default, `renv` only looks for CRAN 
## packages unless we manually specify `bioc::`, but this leads to 
## errors when a package depends on a Bioc package. 
options(repos = BiocManager::repositories())

## Install packages for Rmarkdown compilation
renv::install(c("knitr", "rmarkdown", "BiocStyle"))

## Install SCP.replication from Github
renv::install("remotes")
renv::install("UCLouvain-CBIO/SCP.replication")
## Install latest version of scpdata from Github
renv::install("UCLouvain-CBIO/scpdata")

## Install other packages
renv::install(c("tidyverse", "patchwork", "scuttle", "reticulate", 
                "zellkonverter", "scater", "theislab/destiny", 
                "scran", "ComplexHeatmap", "RColorBrewer", "viridis",
                "circlize", "uwot", "Cairo", "magick", "benchmarkme"))

####---- Load packages ----####

## `renv`’s dependency discovery machinery relies on static analysis 
## of the R code, so it will recognize the following packages using
## the `library` calls. 

library("knitr")
library("rmarkdown")
library("BiocStyle")
library("SCP.replication")
library("SingleCellExperiment")
library("QFeatures")
library("scpdata")
library("scp")
library("tidyverse")
library("patchwork")
library("scuttle")
library("scater")
library("destiny")
library("scran")
library("ComplexHeatmap")
library("RColorBrewer")
library("viridis")
library("circlize")
library("uwot")
library("Cairo")
library("magick")
library("benchmarkme")

####---- Integrate Python environment ----####

library("reticulate")
library("zellkonverter")
renv::use_python() ## Asks a version of Python, selected ~/miniconda3/bin/python3
## `use_python` has created a new directory `python` in the `renv` 
## directory. 

## Prepare the Python environment:
## 1. Install SCeptre (from GitHub), run:
## `pip install https://github.com/bfurtwa/SCeptre/raw/master/Schoof_et_al/code/sceptre-0.1-py3-none-any.whl -t /renv/python/virtualenvs/renv-python-3.8/lib/python3.8/site-packages/`
## 2. Probably that the installer will complain that the numpy version
##    is not compatible. If so, run `pip install "numpy<1.21"`
## 3. Install the `leigenalg` module (missing from the wheel file):
##    `pip install leigenalg`

## Import sceptre to include the modules in renv. The Python 
## dependencies will be stored in the `requirements.txt`. 
reticulate::py_run_string("import sceptre")


####---- Update lockfiles ----####

# renv::status()
renv::snapshot()
## Done! The renv is created and ready to be used. There are two 
## important file: `renv.lock` (containing the list of R package names
## and versions) and `requirements.txt` (containing the list of Python
## module names and versions) 
## 
## To restore the `renv` snapshot in another system, copy the 
## `inst/renvs/schoof2021` directory to the desired system and run:
## `renv::activate("inst/renvs/schoof2021/")`
## `renv::restore()`

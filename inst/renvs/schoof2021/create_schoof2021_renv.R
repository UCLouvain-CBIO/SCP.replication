####---- Initialize renv ----####

## In renv/schoof2021 directory
## It should be only run once
renv::init(bare = TRUE)

####---- Activate renv ----####

## Run this every time you want to use the renv or to modify it
renv::activate("~/PhD/SCP.replication/inst/renvs/schoof2021/")

####---- Install packages ----####

renv::install("BiocManager")
BiocManager::install(version = "devel")
## Add Bioconductor repos. By default, renv only looks for CRAN 
## packages unless we manually specify `bioc::`, but this leads to 
## errors when a package depends on a Bioc package. 
options(repos = BiocManager::repositories())

## Install packages for Rmarkdown compilation
renv::install(c("knitr", "rmarkdown", "BiocStyle"))

## Install SCP.replication from Github (requires dependencies)
renv::install("remotes")
renv::install("UCLouvain-CBIO/SCP.replication")

## Install other packages
renv::install(c("tidyverse", "patchwork", "scuttle", "reticulate", 
                "zellkonverter", "scater", "theislab/destiny", 
                "scran", "ComplexHeatmap", "RColorBrewer", "viridis",
                "circlize", "uwot", "Cairo", "magick", "benchmarkme"))

####---- Load packages ----####

## We need to call `library(PACKAGE,` in order for the package to be 
## detected by renv and stored in the renv.lockfile.

## Note that renv’s dependency discovery machinery relies on static 
## analysis of your R code, and does not understand all of the 
## different ways in which a package might be used in a project. For
## example, renv will detect the following usages:

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
## directory. You need to `pip install` sceptre. To do this:
## 1. Download the wheel file from the GitHub repo: 
## https://github.com/bfurtwa/SCeptre/raw/master/Schoof_et_al/code/sceptre-0.1-py3-none-any.whl
##    https://github.com/bfurtwa/SCeptre/blob/master/Schoof_et_al/code/sceptre-0.1-py3-none-any.whl
## 2. Run `pip install PATH_TO_WHEEL -t [...]/renv/python/virtualenvs/renv-python-3.8/lib/python3.8/site-packages/`
## 3. Probably that the installer will complain that the numpy version
##    is not compatible. If so, run `pip install "numpy<1.21"`
# pip install https://github.com/bfurtwa/SCeptre/raw/master/Schoof_et_al/code/sceptre-0.1-py3-none-any.whl -t /renv/python/virtualenvs/renv-python-3.8/lib/python3.8/site-packages/
    

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
## To do restore this environment snapshot in another system, copy the
## `inst/renvs/schoof2021` directory to the desired system and run:
## `renv::activate("inst/renvs/schoof2021/")`
## `renv::restore()` 

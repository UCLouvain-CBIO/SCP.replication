####---- Initialize renv ----####

wd <- "~/PhD/SCP.replication/inst/renvs/SCoPE2/"
setwd(wd)

## In inst/renvs/SCoPE2 directory
## It should be only run once
renv::init(bare = TRUE)

####---- Activate renv ----####

## Run this every time you want to use the renv or to modify it
renv::activate()

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
renv::install("sva")
renv::install("UCLouvain-CBIO/SCP.replication")

## Install other packages
renv::install(c("tidyverse", "patchwork", "scater", "uwot", 
                "benchmarkme"))

####---- Load packages ----####

## `renv`’s dependency discovery machinery relies on static analysis 
## of the R code, so it will recognize the following packages using
## the `library` calls. 

library("knitr")
library("rmarkdown")
library("BiocStyle")
library("SCP.replication")
library("sva")
library("SingleCellExperiment")
library("QFeatures")
library("scpdata")
library("scp")
library("tidyverse")
library("patchwork")
library("scater")
library("uwot")
library("benchmarkme")

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

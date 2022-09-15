# Single cell replication package

This package contains a series of vignettes that reproduce published
SCP data analyses. Emphasis is put on standardization of the workflows.

## Prerequisites

All vignettes are pre-compiled and can be found on the
[website](https://uclouvain-cbio.github.io/SCP.replication/index.html)
under the `Articles` tab.

The replication package can be installed with R >= 4.1 by running:

```
BiocManager::install("UCLouvain/SCP.replication")
```

## Reproduction of SCP data analyses

The currently available reproduction vignettes are listed below.

#### [Reproduction of the SCoPE2 analysis (Specht et al. 2021)](https://uclouvain-cbio.github.io/SCP.replication/articles/SCoPE2.html)

Project tag: `SCoPE2`

> Specht H, Emmott E, Petelski AA, Huffman RG, Perlman DH, Serra M, et
> al. Single-cell proteomic and transcriptomic analysis of macrophage
> heterogeneity using SCoPE2. Genome Biol. [2021;22:
> 50](http://dx.doi.org/10.1186/s13059-021-02267-5).

See also

- The SCoPE2 [website](https://scope2.slavovlab.net/)
- The SCoPE2 Github [repository](https://github.com/SlavovLab/SCoPE2)

#### [Exploring the autoPOTS data (Liang et al. 2020)](https://uclouvain-cbio.github.io/SCP.replication/articles/liang2020.html)

Project tag: `liang2020`

> Liang, Yiran, Hayden Acor, Michaela A McCown, Andikan J Nwosu,
> Hannah Boekweg, Nathaniel B Axtell, Thy Truong, Yongzheng Cong,
> Samuel H Payne, and Ryan T Kelly. 2020. “Fully Automated Sample
> Processing and Analysis Workflow for Low-Input Proteome Profiling.”
> Anal. Chem., December.2021, [93, 3,
> 1658–1666](https://pubs.acs.org/doi/10.1021/acs.analchem.0c04240).

#### [Reproduction of the hair-cell development analysis (Zhu et al. 2019, eLife)](https://uclouvain-cbio.github.io/SCP.replication/articles/zhu2019EL.html)

Project tag: `zhu2019EL`

> Zhu, Ying, Mirko Scheibinger, Daniel Christian Ellwanger, Jocelyn F
> Krey, Dongseok Choi, Ryan T Kelly, Stefan Heller, and Peter G
> Barr-Gillespie. 2019. “Single-Cell Proteomics Reveals Changes in
> Expression During Hair-Cell Development.” [Elife 8
> (November)](https://elifesciences.org/articles/50777).

#### [Reproduction of the AML model analysis (Schoof et al. 2021)](https://uclouvain-cbio.github.io/SCP.replication/articles/schoof2021.html)

Project tag: `schoof2021`

> Schoof, Erwin M., Benjamin Furtwängler, Nil Üresin, Nicolas Rapin, 
Simonas Savickas, Coline Gentil, Eric Lechman, Ulrich auf Dem Keller, 
John E. Dick, and Bo T. Porse. 2021. “Quantitative Single-Cell 
Proteomics as a Tool to Characterize Cellular Hierarchies.” Nature 
Communications 12 (1): 745679.

See also

- The SCeptre Github [repository](https://github.com/bfurtwa/SCeptre)

#### [Reproducing the multiplexed SCP analysis by Williams et al. 2020](https://uclouvain-cbio.github.io/SCP.replication/articles/williams2020_tmt.html)

Project tag: `williams2020_tmt`

> Williams, Sarah M., Andrey V. Liyu, Chia-Feng Tsai, Ronald J. Moore, 
Daniel J. Orton, William B. Chrisler, Matthew J. Gaffrey, et al. 2020.
“Automated Coupling of Nanodroplet Sample Preparation with Liquid 
Chromatography-Mass Spectrometry for High-Throughput Single-Cell
Proteomics.” Analytical Chemistry 92 (15): 10588–96.

## Replicate the analyses locally

You can reproduce the analysis vignettes on your local machine.
First, clone this repository and `cd` into it. The Rmarkdown files 
that created the articles are available in the `vignettes` directory. 
Then, you need to activate and restore the corresponding `renv` 
project (see tags above). To do so, run the following command: 

```r
install.packages("renv")
renv::activate("inst/renvs/PROJECT_TAG")
renv::restore()
```

Make sure that `renv::restore` completes successfully. This procedure 
installs all the required packages (in R) and modules (in Python) that
are required to run the vignettes. **Warning**: `renv` makes sure to
restore the correct dependency versions, but it does not consider the
R or Python software version. You may run into trouble if your R 
and/or Python versions do not match with the versions indicated at the
end of each analysis (section *Requirements*). 

# Licence

All vignettes are distributed under a 
[CC BY-SA licence](https://creativecommons.org/licenses/by-sa/2.0/) 
licence.
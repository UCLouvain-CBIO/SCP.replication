library(scpdata)
scp <- zhu2019EL()
scp <- filterFeatures(scp,
                      ~ Reverse != "+" &
                          Potential.contaminant != "+")
scp <- zeroIsNA(scp, i = names(scp))
sel <- rowSums(!is.na(assay(scp[["proteins_iBAQ"]]))) >= 5
sel <- names(sel)[sel]

library(biomaRt)
gallus <- useMart("ensembl",
                  dataset = "ggallus_gene_ensembl",
                  host = "https://feb2021.archive.ensembl.org")
ggalusGeneDf <- getBM(values = sel,
                      attributes = c("ensembl_peptide_id", "uniprot_gn_symbol"),
                      filters = "ensembl_peptide_id",
                      mart = gallus)
## Remove duplicate gene symbol
ggalusGeneDf <- ggalusGeneDf[ggalusGeneDf$uniprot_gn_symbol != "CRABP-I", ]
ggalusGeneDf <- ggalusGeneDf[!duplicated(ggalusGeneDf$ensembl_peptide_id), ]

save(ggalusGeneDf, file = "data/ggalusGeneDf.rda")

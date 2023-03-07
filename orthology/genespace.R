#!/usr/bin/env Rscript
library(GENESPACE)
###############################################
# -- change paths to those valid on your system
wd <- "/bigdata/stajichlab/shared/projects/Candida/Candida_lusitaniae_Nanopore/orthology"
path2mcscanx <- "/opt/linux/rocky/8.x/x86_64/pkgs/MCScanX/r51_g97e74f4"
###############################################

gpar <- init_genespace(
  genomeIDs = c("Clus_L17","Clus_ATCC42720"),
  outgroup = NULL,
  ploidy = rep(1,2),
  diamondUltraSens = TRUE,
  wd = wd,
  orthofinderInBlk = FALSE,
  nCores = 4,
  path2orthofinder = "orthofinder",
  path2diamond = "diamond",
  path2mcscanx = path2mcscanx,
)

# -- accomplish the run
out <- run_genespace(gpar,overwrite = T)

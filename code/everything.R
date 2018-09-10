## regenerates all figures from ms/figs/
## runs simulation (sim.R, sim_dep.R)
## and analysis of sdss data (sdss.R,sdss_results.R)
rm(list=ls())

## get rid of old plots/tables
unlink("figs",recursive=TRUE)
unlink("../ms/figs",recursive=TRUE)
dir.create("figs")
dir.create("../ms/figs")

## run simulations and real data
rm(list=ls())
source('sim.R')
rm(list=ls())
source('sim_dep.R')
rm(list=ls())
source('sdss.R')
rm(list=ls())
source('sdss_results.R')

## copy output to ../ms/figs/
rm(list=ls())
fs <- list.files("figs",full.names=TRUE)
file.copy(fs,to="../ms/figs/")

require(mappolymp)
#detach("package:mappoly", unload=TRUE)
require(mappoly)
require(tidyverse)
require(reshape2)
require(plotly)
require(RColorBrewer)
source("~/repos/current_work/mappolymp/R/merge_mappoly_maps.R", echo=TRUE)
map_ch1_BExMG <- readRDS("~/repos/current_work/rose/fullsib_maps/BExMG/map_err_ch_1.rds")
map_ch1_SWxBE <- readRDS("~/repos/current_work/rose/fullsib_maps/SWxBE/map_err_ch_1.rds")
dat_BExMG <- readRDS("~/repos/current_work/rose/data/dat_BExMG.rds")
dat_SWxBE <- readRDS("~/repos/current_work/rose/data/dat_SWxBE.rds")

data.list <- list(BExMG = dat_BExMG,
                  SWxBE = dat_SWxBE)
map.list <- list(BExMG = map_ch1_BExMG,
                 SWxBE = map_ch1_SWxBE)
plot_map_list(map.list)
compare_single_maps(map.list)

parents.mat <- matrix(c("BE", "MG", "SW", "BE"),
                      2, 2, byrow = TRUE,
                      dimnames = list(c("pop1", "pop2"), c("P1", "P2")))

system.time(states <- prepare_maps_to_merge(map.list,
                                            data.list,
                                            parents.mat))

restemp <- hmm_map_reconstruction(ploidy1 = states$ploidy$pl1,
                                  ploidy2 = states$ploidy$pl2,
                                  n.mrk = states$n.mrk,
                                  n.ind = states$n.ind,
                                  haplo = states$states,
                                  emit = states$emit,
                                  rf_vec = rep(0.01, states$n.mrk-1),
                                  verbose = TRUE,
                                  use_H0 = FALSE,
                                  tol = 1e-3)
restemp[[2]]

dummy.map <- map_ch1_BExMG
dummy.map$maps[[1]]$seq.rf <- restemp[[2]]
dummy.map$info$mrk.names <- states$genome.pos$genome.pos$mrk
dummy.map$info$genome.pos <- states$genome.pos$genome.pos$geno.pos

l <- list(SWxBExMG = dummy.map,
          BExMG = map_ch1_BExMG,
        SWxBE = map_ch1_SWxBE)
compare_single_maps(l)

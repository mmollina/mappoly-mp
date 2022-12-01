require(mappoly2)
#detach("package:mappoly", unload=TRUE)
require(mappoly)
require(tidyverse)
require(reshape2)
require(plotly)
source("~/repos/current_work/mappoly2/R/merge_mappoly_maps.R", echo=TRUE)
map_ch1_BExMG <- readRDS("~/repos/current_work/rose/fullsib_maps/BExMG/map_err_ch_1.rds")
map_ch1_SWxBE <- readRDS("~/repos/current_work/rose/fullsib_maps/SWxBE/map_err_ch_1.rds")
dat_BExMG <- readRDS("~/repos/current_work/rose/data/dat_BExMG.rds")
dat_SWxBE <- readRDS("~/repos/current_work/rose/data/dat_SWxBE.rds")

# bla<-ph_matrix_to_list(ph_list_to_matrix(map_ch1_SWxBE$maps[[1]]$seq.ph$Q, 4)[,c(3,1,2,4)])
# names(bla) <- names(map_ch1_SWxBE$maps[[1]]$seq.ph$Q)
# compare_haplotypes(4, bla, map_ch1_SWxBE$maps[[1]]$seq.ph$Q)
# map_ch1_SWxBE$maps[[1]]$seq.ph$Q <- bla

# i <- intersect(map_ch1_BExMG$info$mrk.names, map_ch1_SWxBE$info$mrk.names)
# i1 <- match(i, map_ch1_BExMG$info$mrk.names)[1:200]
# i2 <- match(i, map_ch1_SWxBE$info$mrk.names)[1:200]

i1<-1:500
i2<-1:500

#map_ch1_BExMG <- readRDS("~/repos/current_work/rose/fullsib_maps/BExMG/map_err_ch_1.rds")
#map_ch1_SWxBE <- readRDS("~/repos/current_work/rose/fullsib_maps/SWxBE/map_err_ch_1.rds")
dat <- dat_BExMG
a1 <- get_submap(map_ch1_BExMG, mrk.pos = i1, reestimate.rf = T)
map_ch1_BExMG <- get_submap(map_ch1_BExMG, mrk.pos = i1, reestimate.rf = F)
glimpse(map_ch1_BExMG)
plot(map_ch1_BExMG)

dat <- dat_SWxBE
a2 <- get_submap(map_ch1_SWxBE, mrk.pos = i2, reestimate.rf = T)
map_ch1_SWxBE <- get_submap(map_ch1_SWxBE, mrk.pos = i2, reestimate.rf = F)
#plot(map_ch1_SWxBE)
glimpse(map_ch1_SWxBE)

map_ch1_BExMG$maps[[1]]$seq.rf <- a1$maps[[1]]$seq.rf
map_ch1_SWxBE$maps[[1]]$seq.rf <- a2$maps[[1]]$seq.rf

data.list <- list(BExMG = dat_BExMG,
                  SWxBE = dat_SWxBE)
map.list <- list(BExMG = map_ch1_BExMG,
                 SWxBE = map_ch1_SWxBE)
plot_map_list(map.list)
#compare_single_maps(map.list)

parents.mat <- matrix(c("BE", "MG", "SW", "BE"),
                      2, 2, byrow = TRUE,
                      dimnames = list(c("pop1", "pop2"), c("P1", "P2")))

# system.time(states1 <- prepare_maps_to_merge2(map.list,
#                                  data.list,
#                                  parents.mat))
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
                                  verbose = FALSE,
                                  use_H0 = FALSE,
                                  tol = 1e-3)
dummy.map <- map_ch1_BExMG
dummy.map$maps[[1]]$seq.rf <- restemp[[2]]
dummy.map$info$mrk.names <- states$genome.pos$genome.pos$mrk
dummy.map$info$genome.pos <- states$genome.pos$genome.pos$geno.pos

l <- list(BExMG = map_ch1_BExMG,
        SWxBE = map_ch1_SWxBE,
        SWxBExMG = dummy.map)
compare_single_maps(l)

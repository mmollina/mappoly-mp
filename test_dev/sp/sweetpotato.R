require(mappoly)
require(mappolymp)
require(tidyverse)
require(reshape2)
require(plotly)
source("~/repos/current_work/mappolymp/R/merge_mappoly_maps.R", echo=TRUE)

load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/BT_ch1.rda")
load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/BT_data_ch1.rda")
MAPs_bt <- readRDS("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/MAPs_bt.rds")
map_bt_ch1 <- MAPs_bt[[1]]
plot(map_bt_ch1)

load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/TB_ch1.rda")
load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/TB_data_ch1.rda")
MAPs_tb <- readRDS("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/MAPs_tb.rds")
map_tb_ch1 <- MAPs_tb[[1]]
plot(map_tb_ch1)


load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/NKB_ch1.rda")
load("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/data_by_ch/NKB_data_ch1.rda")
MAPs_nkb <- readRDS("~/repos/my_repos/Beauregard_Tanzania_NewKawogo/BT_TB_NKB_map/full_sib_maps/MAPs_nkb.rds")
map_nkb_ch1 <- MAPs_nkb[[1]]
plot(map_nkb_ch1)


i1<-1:600
i2<-1:360
i3<-1:210

a1 <- get_submap(map_bt_ch1, mrk.pos = i1, reestimate.rf = FALSE)
plot(a1)

a2 <- get_submap(map_tb_ch1, mrk.pos = i2, reestimate.rf = FALSE)
plot(a2)

a3 <- get_submap(map_nkb_ch1, mrk.pos = i3, reestimate.rf = FALSE)
plot(a3)


a1 <- map_bt_ch1
a2 <- map_tb_ch1
a3 <- map_nkb_ch1

map.list <- list(BT = a1, TB = a2, NKB = a3)
compare_single_maps(map.list)

## Beauregard is parent 1 in all pops
id <- Reduce(intersect, lapply(map.list, function(x) x$info$mrk.names))
cor(cbind(BT$dosage.p1[id], TB$dosage.p1[id], NKB$dosage.p1[id]))

parents.mat <- matrix(c("B", "T",
                        "B", "T",
                        "B", "NK"),
                      3, 2, byrow = TRUE,
                      dimnames = list(c("pop1", "pop2", "pop3"), c("P1", "P2")))
data.list <- list(BT, TB, NKB)

dev.off()
system.time(states <- prepare_maps_to_merge(map.list,
                                            data.list,
                                            parents.mat))




x <- mean(unlist(sapply(map.list, function(x) x$maps[[1]]$seq.rf)))
restemp <- hmm_map_reconstruction(ploidy1 = states$ploidy$pl1,
                                  ploidy2 = states$ploidy$pl2,
                                  n.mrk = states$n.mrk,
                                  n.ind = states$n.ind,
                                  haplo = states$states,
                                  emit = states$emit,
                                  rf_vec = rep(x, states$n.mrk-1),
                                  verbose = TRUE,
                                  use_H0 = FALSE,
                                  tol = 10e-3)
restemp[[2]]
dummy.map <- map_bt_ch1
dummy.map$maps[[1]]$seq.rf <- restemp[[2]]/7
dummy.map$info$mrk.names <- states$genome.pos$genome.pos$mrk
dummy.map$info$genome.pos <- states$genome.pos$genome.pos$geno.pos

l <- list(bt_tb_nkb = dummy.map,
          bt = a1,
          tb = a2,
          nkb = a3)
compare_single_maps(l)


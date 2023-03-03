## Gathering parent's phases
w <- table(as.vector(parents.mat))
phases <- vector("list", length(w))
names(phases) <- names(w)

## Gathering ploidy level
pl <- data.list[[1]]$ploidy

## For each unique parent
for(i in names(phases)){
  par.ord <- which(parents.mat == i, arr.ind = T)
  hom.res <- match_homologs(map.list, par.ord, pl)
  L <- ph_matrix_to_list(hom.res$ph)
  names(L) <- rownames(hom.res$ph)
  for(j in 1:nrow(par.ord)){
    ph <- map.list[[par.ord[j,1]]]$maps[[1]]$seq.ph[[par.ord[j,2]]]
    ph.nm <- names(ph)
    ph
  }



  phases[[i]] <- apply(hom.res$ph, 1, paste, collapse = "")
}




bla<-ph_matrix_to_list(ph_list_to_matrix(map_ch1_SWxBE$maps[[1]]$seq.ph$Q, 4)[,c(3,1,2,4)])
names(bla) <- names(map_ch1_SWxBE$maps[[1]]$seq.ph$Q)
compare_haplotypes(4, bla, map_ch1_SWxBE$maps[[1]]$seq.ph$Q)
map_ch1_SWxBE$maps[[1]]$seq.ph$Q <- bla

#i <- intersect(map_ch1_BExMG$info$mrk.names, map_ch1_SWxBE$info$mrk.names)
#i1 <- match(i, map_ch1_BExMG$info$mrk.names)[1:200]
#i2 <- match(i, map_ch1_SWxBE$info$mrk.names)[1:200]

i1<-1:200
i2<-1:200

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




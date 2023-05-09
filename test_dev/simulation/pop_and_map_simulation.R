n.mrk <- 500
ploidy <- 4
h1 <- sim_homologous(ploidy = ploidy,
                     n.mrk = n.mrk,
                     max.d = 2,
                     max.ph = 2,
                     seed = 3467)

h2 <- sim_homologous(ploidy = ploidy,
                     n.mrk = n.mrk,
                     max.d = 2,
                     max.ph = 2,
                     seed = 2379)

P <- list(P1 = h1[c(1,3)],
          P2 = h1[c(2,4)],
          P3 = h2[c(1,3)],
          P4 = h2[c(2,4)])


parents.mat <- matrix(c("P1", "P2",
                        "P2", "P1",
                        "P1", "P3",
                        "P4", "P2",
                        "P2", "P2"),
                      5, 2, byrow = TRUE,
                      dimnames = list(paste0("pop", 1:5),
                                      c("P1", "P2")))

data.list <- vector("list", nrow(parents.mat))
map.list <- vector("list", nrow(parents.mat))
names(data.list) <- names(map.list) <- apply(parents.mat, 1, paste, collapse = "x")
for(i in 1:nrow(parents.mat)){
  Ptemp <- c(P[[parents.mat[i,1]]], P[[parents.mat[i,2]]])
  names(Ptemp) <-  c("hom.allele.p", "p", "hom.allele.q", "q")
  id <- Ptemp$p + Ptemp$q != 0
  Ptemp <- list(hom.allele.p = Ptemp$hom.allele.p[id],
                hom.allele.q = Ptemp$hom.allele.q[id],
                p = Ptemp$p[id],
                q = Ptemp$q[id])
  dat <- poly_cross_simulate(ploidy = ploidy,
                             rf.vec = rep(mf_h(120/n.mrk), n.mrk-1),
                             n.mrk = sum(id),
                             n.ind = 150,
                             hom.allele = Ptemp,
                             seed = i + 2)
  dat$genome.pos <- cumsum(imf_h(c(0,rep(mf_h(120/n.mrk), n.mrk-1))))
  dat$chrom <- rep(1, dat$n.mrk)
  names(dat$genome.pos) <- dat$mrk.names
  names(dat$chrom) <- dat$mrk.names
  dat <- sample_data(input.data = dat, percentage = 70, type = "marker")
  s <- make_seq_mappoly(dat, "all")
  tpt <- est_pairwise_rf(s, ncpus = 30)
  map <- est_rf_hmm_sequential(input.seq = s,
                               start.set = 3,
                               thres.twopt = 10,
                               thres.hmm = 20,
                               extend.tail = 20,
                               info.tail = TRUE,
                               twopt = tpt,
                               sub.map.size.diff.limit = 2,
                               phase.number.limit = 20,
                               reestimate.single.ph.configuration = TRUE,
                               tol = 10e-2,
                               tol.final = 10e-3)
  data.list[[i]] <- dat
  map.list[[i]] <- map
}
save(data.list, map.list, parents.mat,
     file = "~/repos/current_work/mappolymp/test_dev/simulation/data_and_maps.rda")
plot_map_list(map.list, col = mp_pallet1(5))
compare_single_maps(map.list)

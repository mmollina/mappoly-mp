require(mappoly)
h <- list(hom.allele.p = list(c(1), c(1), c(1,2), c(1,2)),
          hom.allele.q = list(c(1), c(1), c(1,2), c(1,2)),
          p = c(1,1,2,2),
          q = c(1,1,2,2))

#h <- sim_homologous(ploidy = 4, n.mrk = 4, max.d = 2, max.ph = 2)

d <- poly_cross_simulate(ploidy = 4,
                         rf.vec = c(0,0.05,0),
                         n.mrk = 4,
                         n.ind = 500,
                         h)
plot(d)
s <- make_seq_mappoly(d, "all")
tpt <- est_pairwise_rf(s)
m <- est_rf_hmm(s,twopt = tpt, tol = 10e-4)
plot(m)

m1 <- get_submap(m, 1:2)
m2 <- get_submap(m, 3:4)
mm <- merge_maps(map.list = list(m1, m2), twopt = tpt, tol = 10e-4)

plot_map_list(list(m,mm))


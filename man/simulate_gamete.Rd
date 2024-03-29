% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiallelic_multiparental_simulation.R
\name{simulate_gamete}
\alias{simulate_gamete}
\title{Gamete simulation}
\usage{
simulate_gamete(homology.group, n.ind, cm.map, prob = NULL, seed = NULL)
}
\arguments{
\item{void}{internal function to be documented}
}
\description{
Gamete simulation
}
\examples{
 ploidy <- 4
 n.mrk <- 7
 n.ind <- 500
 cm.map <- seq(0, 10, length.out = n.mrk)

 h <- simulate_multiallelic_homology_group(ploidy = ploidy,
                                     n.mrk = n.mrk,
                                     alleles = 0:1,
                                     lambda = 1,
                                     shuffle.homolog = TRUE)
 image(h, axes = FALSE, xlab = "Markers", ylab = "homologs")
 abline(h = seq(-1/(2*ploidy - 2), 1 + 1/(2*ploidy - 2),
                length.out = ploidy+1),
        lwd = 2)
 axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h))
 axis(2, at = seq(0,1, length.out = ploidy), labels = colnames(h))
 g <- simulate_gamete(homology.group = h,
                      n.ind = n.ind,
                      cm.map = cm.map)

 dose.geno <- apply(g$gamete, c(1,3), sum)
 colnames(dose.geno) <- paste0("Ind", 1:ncol(dose.geno))
 dat <- data.frame(snp.name = paste0("M", 1:n.mrk),
                   P1 = apply(g$homology.group, 1, sum),
                   P2 = rep(0, n.mrk),
                   sequence = 1,
                   sequence.pos = 1:n.mrk,
                   dose.geno)
 require(mappoly)
 dat <- table_to_mappoly(dat, 4)
 s <- make_seq_mappoly(dat, "all")
 tpt <- est_pairwise_rf(s)
 map <- est_rf_hmm_sequential(s,tpt, tol.final = 10e-5)
 print(map, detailed = TRUE)
 plot(map)
 map
 plot_compare_haplotypes(4,
                         ph_matrix_to_list(g$homology.group[map$info$mrk.names,]),
                         ph_matrix_to_list(g$homology.group[map$info$mrk.names, ]),
                         map$maps[[1]]$seq.ph$P,
                         map$maps[[1]]$seq.ph$Q)

}
\author{
Marcelo Mollinari, \email{mmollin@ncsu.edu}
}

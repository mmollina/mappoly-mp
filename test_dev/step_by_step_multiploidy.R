############## Multiallelism + Multiploidy ################
require(mappoly2)
require(mappoly)
require(tidyverse)
ploidy1 <- 6
ploidy2 <- 2
n.mrk <- 10
n.ind <- 300
al1 <- 1:2
al2 <- 1:2
cm.map <- seq(0, 5, length.out = n.mrk)
####Simulation ####
h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
                                           n.mrk = n.mrk,
                                           alleles = al1,
                                           shuffle.homolog = TRUE)
h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
                                           n.mrk = n.mrk,
                                           alleles = al2,
                                           shuffle.homolog = TRUE)
op <- par(mfrow = c(1,2))
col1 <- sort(unique(as.vector(h1)))
col2 <- sort(unique(as.vector(h2)))
pal <- RColorBrewer::brewer.pal(length(unique(c(col1, col2))), "Set1")
col1 <- pal[col1]
col2 <- pal[col2]
fields::image.plot(h1, axes = FALSE, xlab = "Markers", ylab = "homologs", col = col1)
abline(h = seq(-1/(2*ploidy1 - 2), 1 + 1/(2*ploidy1- 2),
               length.out = ploidy1+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h1))
axis(2, at = seq(0,1, length.out = ploidy1), labels = colnames(h1))
fields::image.plot(h2, axes = FALSE, xlab = "Markers", ylab = "homologs", col = col2)
abline(h = seq(-1/(2*ploidy2 - 2), 1 + 1/(2*ploidy2 - 2),
               length.out = ploidy2+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h2))
axis(2, at = seq(0,1, length.out = ploidy2), labels = colnames(h2))
par(op)
f1 <- simulate_cross(n.ind, h1, h2, cm.map)

####States to visit####
b1 <- apply(f1$ph1, 1, function(x) length(unique(x)))
b2 <- apply(f1$ph2, 1, function(x) length(unique(x)))
info.mrk <- names(which(!(b1==1 & b2==1)))
print(info.mrk)
ngam1 <- choose(ploidy1, ploidy1/2)
ngam2 <- choose(ploidy2, ploidy2/2)
S <- as.matrix(expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1])
n.ind <- dim(f1$offspring)[3]
emit <- h <- vector("list", length(info.mrk))
names(emit) <- names(h) <- info.mrk
emit.temp <- htemp <- vector("list", n.ind)
for(j in info.mrk){
  A1 <- combn(as.numeric(f1$ph1[j,]), ploidy1/2) ## Ordered vector ---> phased
  A2 <- combn(as.numeric(f1$ph2[j,]), ploidy2/2) ## Ordered vector ---> phased
  w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                 apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
  a1 <- apply(str_split_fixed(w,"_", (ploidy1 + ploidy2)/2),1,sort)
  for(i in 1:n.ind){
    id <- which(apply(a1, 2, function(x) all(x == sort(f1$offspring[j,,i]))))
    htemp[[i]] <- as.matrix(S[id, , drop = FALSE])
    emit.temp[[i]] <- matrix(rep(1/nrow(htemp[[i]]), nrow(htemp[[i]])), ncol = 1)
  }
  h[[j]] <- htemp
  emit[[j]] <- emit.temp
}


#### MAPpoly2 - same C++ code as mappoly####
pl1 <- rep(ploidy1, n.ind)
pl2 <- rep(ploidy2, n.ind)
restemp <- hmm_map_reconstruction(ploidy1 = pl1,
                                  ploidy2 = pl2,
                                  n.mrk = length(h),
                                  n.ind = length(h[[1]]),
                                  haplo = h,
                                  emit = emit,
                                  rf_vec = rep(0.01, length(h)-1),
                                  verbose = FALSE,
                                  use_H0 = FALSE,
                                  tol = 10e-5)
(mp.multiallelic <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2))




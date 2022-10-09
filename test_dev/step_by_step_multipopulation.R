require(mappoly2)
require(mappoly)
require(tidyverse)
ploidy.vec <- c(4, 2, 6, 4) #four parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4")
cross.mat <- matrix(c("P1","P2",
                      "P1","P3",
                      "P2","P3",
                      "P1","P1",
                      "P2","P4",
                      "P3","P4"), ncol = 2, byrow = T)
n.ind <- c(60, 60, 60, 60, 60, 60) #per cross
n.mrk <- c(20,15,15,25) #per parent
alleles <- list(P1 = c(1:3),
                P2 = c(4:5),
                P3 = c(6:11),
                P4 = c(12:15))
map.length <- 10 #in centimorgans
sim.cross <- simulate_multiple_crosses(ploidy.vec,
                                       cross.mat,
                                       n.ind,
                                       n.mrk,
                                       alleles,
                                       map.length)
#ploidy1 = pl1,
#ploidy2 = pl2,
#n.mrk = length(h),
#n.ind = length(h[[1]]),
#haplo = h,
#emit = emit,
#rf_vec = rep(0.01, length(h)-1),

####States to visit####
ploidy1 <- sim.cross$pedigree$pl1
ploidy2 <- sim.cross$pedigree$pl2
pd <- sim.cross$pedigree
pd$id <- apply(pd[,3:4], 1, paste0, collapse = "x")
upd <- unique(pd)[,3:4]
unique.cross <- apply(upd, 1, paste0, collapse = "x")
unique.cross.id <- seq_along(unique.cross)
names(unique.cross.id) <- unique.cross
pd$id <- unique.cross.id[pd$id]
upd <- cbind(upd, unique.cross.id)
mrk.names <- rownames(sim.cross$dat[[1]])

for(k in unique.cross.id){
  ngam1 <- choose(upd$pl1[k], upd$pl1[k]/2)
  ngam2 <- choose(upd$pl2[k], upd$pl2[k]/2)
  S <- as.matrix(expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1])
  ###################################
  ###################################
  ###################################
  ########START HERE!!!!!!!!!########
  ###################################
  ###################################
  ###################################
  for(j in mrk.names){
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




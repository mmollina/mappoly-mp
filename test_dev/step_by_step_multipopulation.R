require(mappoly2)
require(mappoly)
require(tidyverse)
ploidy.vec <- c(4, 2, 4, 4) #four parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4")
cross.mat <- matrix(c("P1","P2",
                      "P1","P3",
                      "P2","P3",
                      "P1","P1",
                      "P2","P4",
                      "P3","P4"), ncol = 2, byrow = T)
n.ind <- c(30, 30, 30, 30, 30, 30) #per cross
n.mrk <- c(100,100,100,100) #per parent
alleles <- list(P1 = c(1:4),
                P2 = c(5:6),
                P3 = c(7:10),
                P4 = c(11:14))
map.length <- 10 #in centimorgans
sim.cross <- simulate_multiple_crosses(ploidy.vec,
                                       cross.mat,
                                       n.ind,
                                       n.mrk,
                                       alleles,
                                       map.length)
sim.cross
sim.cross$joint.info
par(bg = "black")
plot(NA, xlim = c(0,10), ylim = c(0,10))
axis(1, col = "white")
axis(2, col = "white")
abline(0,1, col = "red", lwd = 3)
for(it in 1:100){
  cat("it: ", it, "\n")
  sim.cross <- simulate_multiple_crosses(ploidy.vec,
                                         cross.mat,
                                         n.ind,
                                         n.mrk,
                                         alleles,
                                         map.length)
  ####States to visit####
  pd <- sim.cross$pedigree
  pd$id <- apply(pd[,1:2], 1, paste0, collapse = "x")
  upd <- unique(pd)
  unique.cross <- apply(upd[,1:2], 1, paste0, collapse = "x")
  unique.cross.id <- seq_along(unique.cross)
  names(unique.cross.id) <- unique.cross
  pd$id <- unique.cross.id[pd$id]
  upd <- cbind(upd, unique.cross.id)
  mrk.names <- rownames(sim.cross$offspring[[1]])
  ind.names <- rownames(sim.cross$pedigree)
  h <- vector("list", length(mrk.names))
  names(h) <- mrk.names
  for(i in names(h)){
    h[[i]] <- vector("list", length(ind.names))
    names(h[[i]]) <- ind.names
  }
  emit <- h
  for(k in unique.cross.id){
    ngam1 <- choose(upd$pl1[k], upd$pl1[k]/2)
    ngam2 <- choose(upd$pl2[k], upd$pl2[k]/2)
    S <- as.matrix(expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1])
    P1 <- upd[k,"Par1"]
    P2 <- upd[k,"Par2"]
    pl1 <- upd$pl1[k]
    pl2 <- upd$pl2[k]
    cur.ind.names <- rownames(pd)[pd$id == k]
    for(j in mrk.names){
      a1 <- sim.cross$phases[[P1]][j,]
      a2 <- sim.cross$phases[[P2]][j,]
      if(any(is.na(a1)) | any(is.na(a2))){
        for(i in cur.ind.names){
          h[[j]][[i]] <- S
          emit[[j]][[i]] <- matrix(rep(1/nrow(h[[j]][[i]]), nrow(h[[j]][[i]])), ncol = 1)
        }
      } else {
        A1 <- combn(a1, pl1/2) ## Ordered vector ---> phased
        A2 <- combn(a2, pl2/2) ## Ordered vector ---> phased
        w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                       apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
        u <- str_split_fixed(w,"_", (pl1 + pl2)/2)
        storage.mode(u) <- "integer"
        v <- apply(u,1,sort)
        for(i in cur.ind.names){
          id <- which(apply(v, 2, function(x) all(x == sort(sim.cross$offspring[[i]][j,]))))
          h[[j]][[i]] <- S[id, , drop = FALSE]
          emit[[j]][[i]] <- matrix(rep(1/nrow(h[[j]][[i]]), nrow(h[[j]][[i]])), ncol = 1)
        }
      }
    }
  }
  #### MAPpoly2 - same C++ code as mappoly####
  system.time(restemp <- hmm_map_reconstruction(ploidy1 = sim.cross$pedigree$pl1,
                                    ploidy2 = sim.cross$pedigree$pl2,
                                    n.mrk = length(h),
                                    n.ind = length(h[[1]]),
                                    haplo = h,
                                    emit = emit,
                                    rf_vec = rep(0.01, length(h)-1),
                                    verbose = FALSE,
                                    use_H0 = FALSE,
                                    tol = 1e-3))
  x <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
  y <- round(sim.cross$map$map, 2)
  z <- lm(y~x)
  points(y~x, cex = .5, pch = 20, col = "pink")
  abline(z, col = "yellow", lwd = .5)
}








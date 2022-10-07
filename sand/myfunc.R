myfunc <- function(ploidy1, ploidy2, al1, al2, n.mrk, n.ind, cm.map){
  ####Simulation ####
  h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
                                             n.mrk = n.mrk,
                                             alleles = al1,
                                             shuffle.homolog = TRUE)
  h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
                                             n.mrk = n.mrk,
                                             alleles = al2,
                                             shuffle.homolog = TRUE)

  f1 <- simulate_cross(n.ind, h1, h2, cm.map)

  b1 <- apply(f1$ph1, 1, function(x) length(unique(x)))
  b2 <- apply(f1$ph2, 1, function(x) length(unique(x)))
  info.mrk <- names(which(!(b1==1 & b2==1)))
  ####States to visit####
  ngam1 <- choose(ploidy1, ploidy1/2)
  ngam2 <- choose(ploidy2, ploidy2/2)
  S <- as.matrix(expand.grid(0:(ngam1-1), 0:(ngam2-1))[,2:1])
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
  restemp <- hmm_map_reconstruction(ploidy = ploidy1,
                                    n.mrk = length(h),
                                    n.ind = length(h[[1]]),
                                    haplo = h,
                                    emit = emit,
                                    rf_vec = rep(0.01, length(h)-1),
                                    verbose = FALSE,
                                    use_H0 = FALSE,
                                    tol = 10e-5)
  mpm <- cumsum(mappoly::imf_h(c(0, restemp[[2]])))
  names(mpm) <- info.mrk
  return(round(mpm,2))
}

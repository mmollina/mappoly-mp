myfunc <- function(){
  ####Simulation ####
  h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
                                             n.mrk = n.mrk,
                                             alleles = al1,
                                             shuffle.homolog = TRUE)
  h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
                                             n.mrk = n.mrk,
                                             alleles = al2,
                                             shuffle.homolog = TRUE)

  f1 <- simulate_cross(h1, h2, cm.map)

  ####States to visit####
  ngam1 <- choose(ploidy1, ploidy1/2)
  ngam2 <- choose(ploidy2, ploidy2/2)
  S <- cbind(ploidy1, ploidy2, 1,
             expand.grid(0:(ngam1-1), 0:(ngam2-1))[,2:1],
             0:((ngam1*ngam2)-1), 0)
  colnames(S) <- c("pl1", "pl2", "pl.id", "st.p1", "st.p2", "st.all", "emit")
  Y <- NULL
  b1 <- apply(f1$ph1, 1, function(x) length(unique(x)))
  b2 <- apply(f1$ph2, 1, function(x) length(unique(x)))
  info.mrk <- as.numeric(str_remove_all(names(which(!(b1==1 & b2==1))), "M"))
  for(i in 1:dim(f1$offspring)[3]){
    for(j in info.mrk){
      A1 <- combn(as.numeric(f1$ph1[j,]), ploidy1/2) ## Ordered vector ---> phased
      A2 <- combn(as.numeric(f1$ph2[j,]), ploidy2/2) ## Ordered vector ---> phased
      w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                     apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
      a1 <- apply(str_split_fixed(w,"_", (ploidy1 + ploidy2)/2),1,sort)
      id <- which(apply(a1, 2, function(x) all(x == sort(f1$offspring[j,,i]))))
      tmp <- S[id,]
      tmp$emit[] <- log(1/length(tmp$emit))
      Y<-rbind(Y, cbind(ind = paste0("Ind_", i), mrk = paste0("M_", j), tmp))
    }
  }
  sts <- list(hmm.info = Y,
              err = 0,
              is.log = TRUE,
              ploidy.cross.id = matrix(c(ploidy1,ploidy2),1,2, dimnames = list(paste(ploidy1, ploidy2, sep = "x"), NULL)))
  mrks <- as.character(unique(sts$hmm.info$mrk))
  emit <- h <- vector("list", length(mrks))
  names(emit)  <- names(h) <- mrks
  for(i in mrks)
  {
    st.temp <- sts$hmm.info %>%
      filter(mrk == i)
    ind <- as.character(unique(st.temp$ind))
    emit.temp <- htemp <- vector("list", length(ind))
    names(emit.temp) <- names(htemp) <- ind
    for(j in ind){
      htemp[[j]] <- as.matrix(st.temp %>%
                                filter(ind == j) %>%
                                select("st.p1", "st.p2"))
      emit.temp[[j]] <- exp(as.matrix(st.temp %>%
                                        filter(ind == j) %>%
                                        select("emit")))
    }
    h[[i]] <- htemp
    emit[[i]] <- emit.temp
  }

  #### MAPpoly2 - same C++ code as mappoly####
  restemp <- hmm_map_reconstruction(ploidy = ploidy,
                                    n.mrk = length(h),
                                    n.ind = length(h[[1]]),
                                    haplo = h,
                                    emit = emit,
                                    rf_vec = rep(0.01, length(h)-1),
                                    verbose = FALSE,
                                    use_H0 = FALSE,
                                    tol = 10e-5)
  mp.multiallelic <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
  list(mp.multiallelic, cm.map)
}

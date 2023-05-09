

#' Expected genotypic frequency
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
segreg_multiallele <- function(x,y)
{
  z2 <- z1 <- NULL
  for(i in 1:length(x)){
    z1<-c(z1, rep(names(x)[i], x[i]))
    z2<-c(z2, rep(names(y)[i], y[i]))
  }
  z1<-combn(z1, length(z1)/2)
  z2<-combn(z2, length(z2)/2)
  id<-expand.grid(1:ncol(z1), 1:ncol(z1))
  w <- character(nrow(id))
  for(i in 1:nrow(id))
    w[i] <- paste0(sort(c(z1[,id[i,1]], z2[,id[i,2]])), collapse = "")
  w<-table(w)
  return(w/sum(w))
}
#' Short to long format
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
short_to_long<-function(x)
{
  if(any(is.na(x))) return(NA)
  paste0(sort(rep(names(x), x)), collapse = "")
}
#' Chisq test for a single SNP
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom XNomial "xmonte"
#' @export
segreg_test_single_snp <- function(dat, mrk.id)
{
  m<- dat$m
  p <- segreg_multiallele(dat$dosage.p[[mrk.id]], dat$dosage.q[[mrk.id]])
  d <- dat$geno.dose[[mrk.id]]
  res<-character(nrow(d))
  for(i in 1:nrow(d)){
    res[i]<-short_to_long(d[i,])
  }
  res[nchar(res) > m] <- NA
  res <- table(res)
  x <- p
  x[] <- 0
  x[names(res)] <- res
  x<-x[names(p)]
  if(sum(x) == 0) return(NA)
  if(length(x) > 1){
    out <- XNomial::xmonte(x, p)
    return(out$pLLR)
  }
  else return(NA)
}
#' Chisq test for the whole population
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
segreg_test_multiallelic <- function(dat)
{
  p.val<-rep(NA, dat$n.mrk)
  for(i in 1:dat$n.mrk){
    #cat(".")
    #if(i%%50==0) cat("\n")
    suppressWarnings(ptemp<-segreg_test_single_snp(dat, i))
    if(is.na(ptemp)) next()
    p.val[i] <- ptemp
  }
  cat("\n")
  names(p.val) <- dat$mrk.names
  list(p.val = p.val, allelic.level = sapply(sapply(dat$dosage.p, names), length))
}
#' Plot chisq test for the whole population
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme element_text scale_color_manual
#' @export
plot_segreg_multiallelic <- function(p.val)
{
  df <- data.frame(mrk.name = names(p.val$p.val), allelic.level = as.factor(p.val$allelic.level),
                   p.val = -log10(p.val$p.val))
  ggplot(df, aes(mrk.name, p.val, colour = allelic.level)) +
    geom_point() + geom_hline(aes(yintercept = -log10(0.01))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) #+
    #scale_color_manual(values=my_pal)
}
#' Returns the id of states that should be visited
#' in a full state HMM
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_states_to_visit <- function(dat, k)
{
  m<-dat$m
  x <- dat$dosage.p[[k]]
  y <- dat$dosage.q[[k]]
  z <- apply(dat$geno.dose[[k]], 1, short_to_long)
  ## one possible permutation
  z1 <- strsplit(short_to_long(x), "")[[1]]
  z2 <- strsplit(short_to_long(y), "")[[1]]
  z1<-combn(z1, length(z1)/2)
  z2<-combn(z2, length(z2)/2)
  id<-expand.grid(1:ncol(z1), 1:ncol(z1))[,2:1]
  w <- character(nrow(id))
  for(i in 1:nrow(id))
    w[i] <- paste0(sort(c(z1[,id[i,1]], z2[,id[i,2]])), collapse = "")
  ngam <- choose(m, m/2)
  A<-as.matrix(expand.grid(0:(ngam-1),
                           0:(ngam-1))[,2:1])
  ## states to visit
  h<-vector("list", length(z))
  for(i in 1:length(h)){
    h[[i]] <- A[which(w %in% z[i]), , drop = FALSE]
  }
  h
}
#' Returns the id of states that should be visited
#' in a full state HMM over all markers in a data
#' set
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
states_for_all_markers <- function(dat)
{
  h<-vector("list", dat$n.mrk)
  for(k in 1:length(h)){
    cat(k, "\n")
    h[[k]] <- get_states_to_visit(dat, k)
  }
  h
}

#' Estimate a genetic map given a sequence of block markers
#'
#' @param void internal function to be documented
#' @export est_haplo_hmm
est_haplo_hmm <-
  function(m, n.mrk, n.ind, haplo, emit = NULL,
           rf_vec, verbose = TRUE, use_H0 = FALSE,
           highprec = FALSE, tol = 10e-4) {
    ## Checking capabilities
    if (verbose && !capabilities("long.double") && highprec){
      cat("This function uses high precision calculations, but your system's architecture doesn't support long double allocation ('capabilities('long.double') = FALSE'). Running in low precision mode.\n")
      highprec = FALSE
    }
    ## In case no genotypic probabilities distributions are provided
    if(is.null(emit)){
      emit <- vector("list", length(haplo))
      for(i in  1:length(haplo)){
        tempemit <- vector("list", length(haplo[[i]]))
        for(j in 1:length(haplo[[i]])){
          tempemit[[j]] <- rep(1, nrow(haplo[[i]][[j]]))
        }
        emit[[i]] <- tempemit
      }
    }
    if(highprec){
     # res.temp <-
    #    .Call("est_haplotype_map_highprec",
    #          m,
    ##          n.mrk,
    #          n.ind,
    #          haplo,
    #          emit,
    ##          rf_vec,
    #          verbose,
    #          tol,
    #          use_H0,
    #          PACKAGE = "mappoly")
    #  return(res.temp)

    } else {
      #res.temp <-
      #  .Call("est_haplotype_map",
      #        m,
      #        n.mrk,
      #        n.ind,
      ##        haplo,
      #        emit,
      #        rf_vec,
      #        verbose,
      #        tol,
      #        use_H0,
      #        PACKAGE = "mappolymp")
      #return(res.temp)
    }
  }


#' Returns the recombination fraction matrix for
#' a multiallelic set of markers
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
multiallelic_rf_mat <- function(dat,h)
{
  M<-matrix(NA, dat$n.mrk, dat$n.mrk)
  for(i in 1:(dat$n.mrk-1)){
    for(j in i:dat$n.mrk){
      cat(i, " ", j, "\n")
      htemp<-h[c(i,j)]
      res <- est_haplo_hmm(m = dat$m, n.mrk = length(htemp), n.ind = dat$n.ind, haplo = htemp,
                                    rf_vec = rep(0.03, length(htemp)-1), verbose = F,
                                    use_H0 = FALSE, tol = 10e-4)
      M[i,j] <- M[j,i] <- res[[2]]
    }
  }
  M
}
#' Given multialleleic marker information, it returns
#' a list of phased sub-maps
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom mappoly ph_matrix_to_list
#' @export
multi_2_hapmap <- function(p,q,g,m, mrk.names)
{
  M1 <- M2 <- matrix(0, length(p), m, dimnames = list(0:(length(p)-1), 1:m))
  j1<-j2<-1
  for(i in names(p)){
    if(p[i] != 0){
      M1[i,j1:(p[i]+j1-1)] <- 1
      j1 <- j1 + p[i]
    }
    if(q[i] != 0){
      M2[i,j2:(q[i]+j2-1)] <- 1
      j2 <- j2 + q[i]
    }
  }
  map <- structure(list(info = list(m = m,
                                    n.mrk = length(mrk.names),
                                    data.name = "",
                                    seq.num  = seq_along(mrk.names),
                                    mrk.names = mrk.names),
                        maps = list(list(seq.ph = list(P = ph_matrix_to_list(M1),
                                                       Q = ph_matrix_to_list(M2)),
                                         seq.rf = rep(1e-7, length(mrk.names)-1),
                                         seq.num = seq_along(mrk.names),
                                         loglike = 0)
                        )
  ), class = "mappoly.map"
  )
  dp <- apply(M1, 1, sum)
  dq <- apply(M2, 1, sum)
  place.holder <- rep(NA, map$info$n.mrk)
  names(place.holder) <- names(dp)
  dattemp <- structure(list(m = m,
                            n.ind = nrow(g),
                            n.mrk = map$info$n.mrk,
                            ind.names = rownames(g),
                            mrk.names = map$info$mrk.names,
                            dosage.p = dp,
                            dosage.q = dq,
                            sequence = place.holder,
                            sequence.pos = place.holder,
                            seq.ref = NULL,
                            seq.alt = NULL,
                            all.mrk.depth = place.holder,
                            prob.thres = NA,
                            geno.dose = t(g[,names(dp)]),
                            nphen = 0,
                            phen = NULL,
                            kept = NULL,
                            elim.correspondence = NULL,
                            chisq.pval = place.holder),
                       class = "mappoly.data")
  rownames(dattemp$geno.dose) <- dattemp$mrk.names
  return(list(map = map, dat = dattemp))
}
#' Multiallelic to biallelic data
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @importFrom mappoly merge_datasets
#' @export
multi_2_biallelic <- function(dat)
{
  m<-dat$m
  haplo <- vector("list", dat$n.mrk)
  biallelic.dat <- NULL
  for(i in 1:dat$n.mrk){
    cat("mrk", i, "\n")
    p <- dat$dosage.p[[i]]
    q <- dat$dosage.q[[i]]
    g <- dat$geno.dose[[i]]
    a<-multi_2_hapmap(p,q,g,m, mrk.names = paste0("M_", i, "_", 1:length(p)))
    a$map$info$data.name <- "biallelic.dat"
    haplo[[i]]<-a$map
    if(i != 1){
      names(biallelic.dat)[c(1,6:9)] <- c("ploidy", "dosage.p1", "dosage.p2", "chrom", "genome.pos")
      names(a$dat)[c(1,6:9)] <- c("ploidy", "dosage.p1", "dosage.p2", "chrom", "genome.pos")
      haplo[[i]]$info$seq.num <- haplo[[i]]$info$seq.num + max(haplo[[i-1]]$info$seq.num)
    }
    haplo[[i]]$maps[[1]]$seq.num <- haplo[[i]]$info$seq.num
    biallelic.dat <- mappoly::merge_datasets(biallelic.dat, a$dat)
  }
  list(biallelic.dat = biallelic.dat, haplotypes = haplo)
}
#' Get states to visit in a pair of markers
#'
#' @param void interfunction to be documented
#' @keywords internal
#' @export
get_states_to_visit_pair <- function(input.map1,
                                     input.map2,
                                     dat,
                                     pos,
                                     rf.matrix)
{
  m<-dat$m
  d<-dat$geno.dose[[pos]]
  z <- apply(d[,sort(colnames(d))], 1, short_to_long)
  y <- generate_all_link_phases_elim_equivalent_haplo(block1 = input.map1$maps[[1]],
                                                      block2 = input.map2$maps[[1]],
                                                      rf.matrix = rf.matrix,
                                                      ploidy = m,
                                                      max.inc = 0)
  ALL.H <- vector("list", length(y))
  for(k in 1:length(y)){
    z1 <- as.character(apply(ph_list_to_matrix(y[[k]]$P, m), 2, function(x) which(x==1))-1)
    z2 <- as.character(apply(ph_list_to_matrix(y[[k]]$Q, m), 2, function(x) which(x==1))-1)
    z1<-combn(z1, length(z1)/2)
    z2<-combn(z2, length(z2)/2)
    id<-expand.grid(1:ncol(z1), 1:ncol(z1))[,2:1]
    w <- character(nrow(id))
    for(i in 1:nrow(id))
      w[i] <- paste0(sort(c(z1[,id[i,1]], z2[,id[i,2]])), collapse = "")
    ngam <- choose(m, m/2)
    A<-as.matrix(expand.grid(0:(ngam-1),
                             0:(ngam-1))[,2:1])
    ## states to visit
    h<-vector("list", length(z))
    for(i in 1:length(h)){
      h[[i]] <- A[which(w %in% z[i]), , drop = FALSE]
    }
    ALL.H[[k]] <- h
  }
  list(states = ALL.H, phases = y)
}

#' Simulate multiallelic homology groups
#'
#' @param void internal function to be documented
#' @examples
#'  ## Multiallelic
#'  h1 <- simulate_multiallelic_homology_group(ploidy = 4,
#'                                     n.mrk = 7,
#'                                     alleles = 3:7,
#'                                     lambda = 10)
#'  ## Biallelic with shuffling
#'  h2 <- simulate_multiallelic_homology_group(ploidy = 4,
#'                                      n.mrk = 7,
#'                                      alleles = 0:1,
#'                                      lambda = 1,
#'                                      shuffle.homolog = TRUE)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
simulate_multiallelic_homology_group <- function(ploidy,
                                          n.mrk,
                                          alleles = 0:(ploidy-1),
                                          lambda = median(alleles),
                                          shuffle.homolog = FALSE,
                                          seed = NULL)
{
  if(!is.null(seed)) set.seed(seed)
  P <- matrix(NA, n.mrk, ploidy)
  for(i in 1:n.mrk)
    P[i,] <- sort(sample(alleles,
                         size = ploidy,
                         prob = dpois(x = alleles, lambda = lambda),
                         replace = TRUE))

  if(shuffle.homolog){
    P <- t(apply(P, 1, sample))
  } else{
    P <- t(apply(P, 1, sort))
  }
  dimnames(P) <- list(paste0("M", 1:nrow(P)), paste0("H", 1:ncol(P)))
  P
}
#' Gamete simulation
#'
#' @param void internal function to be documented
#' @examples
#'  ploidy <- 4
#'  n.mrk <- 7
#'  n.ind <- 500
#'  cm.map <- seq(0, 10, length.out = n.mrk)
#'
#'  h <- simulate_multiallelic_homology_group(ploidy = ploidy,
#'                                      n.mrk = n.mrk,
#'                                      alleles = 0:1,
#'                                      lambda = 1,
#'                                      shuffle.homolog = TRUE)
#'  image(h, axes = FALSE, xlab = "Markers", ylab = "homologs")
#'  abline(h = seq(-1/(2*ploidy - 2), 1 + 1/(2*ploidy - 2),
#'                 length.out = ploidy+1),
#'         lwd = 2)
#'  axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h))
#'  axis(2, at = seq(0,1, length.out = ploidy), labels = colnames(h))
#'  g <- simulate_gamete(homology.group = h,
#'                       n.ind = n.ind,
#'                       cm.map = cm.map)
#'
#'  dose.geno <- apply(g$gamete, c(1,3), sum)
#'  colnames(dose.geno) <- paste0("Ind", 1:ncol(dose.geno))
#'  dat <- data.frame(snp.name = paste0("M", 1:n.mrk),
#'                    P1 = apply(g$homology.group, 1, sum),
#'                    P2 = rep(0, n.mrk),
#'                    sequence = 1,
#'                    sequence.pos = 1:n.mrk,
#'                    dose.geno)
#'  require(mappoly)
#'  dat <- table_to_mappoly(dat, 4)
#'  s <- make_seq_mappoly(dat, "all")
#'  tpt <- est_pairwise_rf(s)
#'  map <- est_rf_hmm_sequential(s,tpt, tol.final = 10e-5)
#'  print(map, detailed = TRUE)
#'  plot(map)
#'  map
#'  plot_compare_haplotypes(4,
#'                          ph_matrix_to_list(g$homology.group[map$info$mrk.names,]),
#'                          ph_matrix_to_list(g$homology.group[map$info$mrk.names, ]),
#'                          map$maps[[1]]$seq.ph$P,
#'                          map$maps[[1]]$seq.ph$Q)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom gtools permutations
#' @export
simulate_gamete <- function(homology.group,
                            n.ind,
                            cm.map,
                            prob = NULL,
                            seed = NULL)
{
  ploidy <- ncol(homology.group)
  n.mrk <- nrow(homology.group)
  if(!is.null(seed)) set.seed(seed)
  dist.vec <- diff(cm.map)
  rf.vec <-   0.5 * (1 - exp(-dist.vec/50)) #Haldane
  if(length(rf.vec) == 1) rf.vec <- rep(rf.vec, n.mrk-1)
  res <- array(NA, c(n.mrk, ploidy/2 ,n.ind))
  rf.res <- numeric(n.mrk-1)

  ## Listing all possible bivalent configurations
  a <- permutations(ploidy, ploidy, 1:ploidy)
  bv.conf <- vector("list", nrow(a))
  for(i in 1:nrow(a))
  {
    temp <- apply(matrix(a[i,], 2, ploidy/2), 2, sort)
    bv.conf[[i]] <- temp[,order(temp[1,]), drop = FALSE]
  }
  bv.conf <- unique(bv.conf)
  if(is.null(prob))
    prob <- rep(1/length(bv.conf), length(bv.conf))

  for(k in 1:n.ind){              #for each individual
    gen.1 <- matrix(1:ploidy,ploidy,n.mrk)  #completely informative markers in 'n.mrk' positions
    choosed_biv <- sample(bv.conf, 1, prob = prob)[[1]] #sampling one bivalent configuration based on given probabilities
    for(i in 1:ncol(choosed_biv))
      choosed_biv[,i] <- sample(choosed_biv[,i])
    pole.1 <- choosed_biv[1,,drop = FALSE]
    pole.2 <- choosed_biv[2,,drop = FALSE]
    set.2 <- gen.1[pole.1,,drop = FALSE]      #allocating the chromosomes on the variables set.1 and set.2, thus (set.1[i], set.2[i]) represents a bivalent
    set.1 <- gen.1[pole.2,,drop = FALSE]
    for(i in 1:(ploidy/2)){         #for each one of the ploidy/2 chromosome pair (bivalents)
      a <- set.1[i,]
      b <- set.2[i,]
      for(j in 1:(n.mrk-1)){             #for each adjacent interval between.markers
        if(runif(1)  < rf.vec[j]){       #if a random number drawn from the interval [0,1] (according a uniform distribution)
          #is less than the recombination fraction for that interval
          which.swap <- c((j+1):n.mrk)     #the alleles for that interval and bivalent are swapped
          temp <- a[which.swap]
          a[which.swap] <- b[which.swap]
          b[which.swap] <- temp
        }
      }             #this completes the whole bivalent
      set.1[i,] <- a  #attributing the resulting vector to the initial variables
      set.2[i,] <- b
    }               #for all bivalents

    if(sample(0:1,1)) {
      gam <- set.1 } else gam <- set.2 #sample one of the meiotic products
    for(i in 1:(ploidy/2)){ #counting the recombinant chromosomes in their multiallelic form
      for(j in 2:ncol(gam)){
        if(!gam[i,j] == gam[i,j-1])
          rf.res[j-1] <- rf.res[j-1]+1
      }
    }
    for(i in 1:n.mrk)
      res[i,,k] <- as.numeric(homology.group[i,gam[,i]])
  }
  rf.calc <- rf.res/(n.ind*ploidy/2)  #computing the recombination fraction
  list(gamete = res, c.o.count = rf.res, homology.group = homology.group)
}

#' Simulate a biparental cross
#'
#' @param void internal function to be documented
#' @examples
#'  ploidy1 <- 4
#'  ploidy2 <- 4
#'  n.mrk <- 7
#'  n.ind <- 500
#'  cm.map <- seq(0, 10, length.out = n.mrk)
#'  h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
#'                                             n.mrk = n.mrk,
#'                                             alleles = 0:1,
#'                                             lambda = 1,
#'                                             shuffle.homolog = F)
#'  op <- par(mfrow = c(1,2))
#'  image(h1, axes = FALSE, xlab = "Markers", ylab = "homologs")
#'  abline(h = seq(-1/(2*ploidy1 - 2), 1 + 1/(2*ploidy1- 2),
#'                 length.out = ploidy1+1),
#'         lwd = 2)
#'  axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h1))
#'  axis(2, at = seq(0,1, length.out = ploidy1), labels = colnames(h1))
#'
#'  h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
#'                                             n.mrk = n.mrk,
#'                                             alleles = 0:1,
#'                                             lambda = 1,
#'                                             shuffle.homolog = F)
#'  image(h2, axes = FALSE, xlab = "Markers", ylab = "homologs")
#'  abline(h = seq(-1/(2*ploidy2 - 2), 1 + 1/(2*ploidy2 - 2),
#'                 length.out = ploidy2+1),
#'         lwd = 2)
#'  axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h2))
#'  axis(2, at = seq(0,1, length.out = ploidy2), labels = colnames(h2))
#'  par(op)
#'  f1 <- simulate_cross(h1, h2, cm.map)
#'
#'  dose.geno <- apply(f1$offspring, c(1,3), sum)
#'  colnames(dose.geno) <- paste0("Ind", 1:ncol(dose.geno))
#'  dat <- data.frame(snp.name = paste0("M", 1:n.mrk),
#'                    P1 = apply(f1$ph1 , 1, sum),
#'                    P2 = apply(f1$ph2 , 1, sum),
#'                    sequence = 1,
#'                    sequence.pos = 1:n.mrk,
#'                    dose.geno)
#'  require(mappoly)
#'  dat <- table_to_mappoly(dat, 4)
#'  s <- make_seq_mappoly(dat, "all")
#'  tpt <- est_pairwise_rf(s)
#'  map <- est_rf_hmm_sequential(s,tpt, tol.final = 10e-5)
#'  print(map, detailed = TRUE)
#'  plot(map)
#'  plot_compare_haplotypes(4,
#'                          ph_matrix_to_list(f1$ph1[map$info$mrk.names,]),
#'                          ph_matrix_to_list(f1$ph2[map$info$mrk.names,]),
#'                          map$maps[[1]]$seq.ph$P,
#'                          map$maps[[1]]$seq.ph$Q)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
simulate_cross <- function(n.ind,
                           h1,
                           h2,
                           cm.map,
                           prob1 = NULL,
                           prob2 = NULL,
                           seed = NULL){
  if(!is.null(seed)) set.seed(seed)
  ##Parent 1
  data.P1 <- simulate_gamete(h1, n.ind, cm.map, prob1)
  ##Parent 2
  data.P2  <- simulate_gamete(h2, n.ind, cm.map, prob2)
  rf.calc <- (data.P1$c.o.count + data.P2$c.o.count)/(n.ind*(ncol(h1)+ncol(h2))/2)
  offspring <- array(NA, dim = c(nrow(h1), (ncol(h1)+ncol(h2))/2, n.ind),
                     dimnames = list(rownames(h1), NULL, paste0("Ind_", 1:n.ind)))
  for(i in 1:n.ind)
    offspring[,,i] <- cbind(data.P1$gamete[,,i], data.P2$gamete[,,i])
  list(offspring = offspring,
       rf.calc = rf.calc,
       ph1 = data.P1$homology.group,
       ph2 = data.P2$homology.group)
}

#' Simulate haplotypes in connected populations
#'
#' @param void internal function to be documented
#' @examples
#' h.c <- simulate_connected(ploidy = c(4,2,6),
#'                            n.mrk = c(10, 15, 10),
#'                            alleles = list(c(0:3), c(3:4), c(3:6)),
#'                            map.length = 10)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom stringr str_split str_split_fixed
#' @importFrom tibble rownames_to_column
#' @importFrom plyr join_all
#' @export
simulate_connected <- function(ploidy,
                                n.mrk,
                                alleles,
                                map.length,
                                lambda = sapply(alleles, median)){
  shared <- NULL
  for(i in length(ploidy):1){
    v <- combn(1:length(ploidy), i)
    shared <- c(shared,apply(v, 2, paste0, collapse = ""))
  }
  v <- str_split(shared, "")
  v <- v[sapply(v, length) > 1]
  p <- table(unlist(v))
  h <- vector("list", length(ploidy))
  for(i in 1:length(v)){
    a <- min(ceiling(n.mrk/p))
    map <- sort(runif(a, 0, map.length))
    names(map) <- paste0(shared[i],"_", 1:length(map))
    for(j in 1:length(ploidy)){
      if(stringr::str_detect(shared[i], as.character(j))){
        h.t <- simulate_multiallelic_homology_group(ploidy[j], length(map), alleles[[j]], lambda[j])
        if(j == 1)
           h.t <- matrix(c(1,1,1,1,1,0,0,0,0,0,0,0), ncol = 4)
        else
           h.t[] <- 0
        colnames(h.t) <- paste0("P",j,"_", 1:ploidy[j])
        h.t <- cbind(h.t, map)
        h[[j]] <- rbind(h[[j]], h.t)
      }
    }
  }
  res <- lapply(h, function(x) x[order(x[,ncol(x)]),])
  res <- lapply(res, function(x) rownames_to_column(as.data.frame(x), "mrks"))
  res <- join_all(res, by='map', type = 'full')
  res <- res[,sort(colnames(res))]
  res <- res[order(res$map),]
  res
}

#' Simulate multiple crosses based
#'
#' @param void internal function to be documented
#' @examples
#' p <- c(4, 2, 6, 4)
#' names(p) <- c("P1", "P2", "P3", "P4")
#' cm <- matrix(c("P1","P2","P1","P3","P2","P3","P1","P1","P2","P4","P3","P4"), ncol = 2, byrow = T)
#' D <- simulate_multiple_crosses(ploidy = p, #four parents
#'                                cross.mat = cm,
#'                                n.ind = c(60, 60, 60, 60, 60, 60), #per cross
#'                                n.mrk= c(100,100,100,100), # per parent
#'                                alleles = list(c(0:3), c(4:5), c(6:11), c(12:15)),
#'                                map.length = 100)
#' D
#' names(D)
#' names(D$dat)
#' D$dat$`1x2`
#' dimnames(D$dat$`1x2`)
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom stringr str_detect str_remove_all
#' @export
simulate_multiple_crosses <- function(ploidy,
                                      cross.mat,
                                      n.ind,
                                      n.mrk,
                                      alleles,
                                      map.length,
                                      lambda = sapply(alleles, median))
{
  h <- simulate_connected(ploidy = ploidy,
                           n.mrk = n.mrk,
                           map.length = map.length,
                           alleles = alleles,
                           lambda = lambda)
  rownames(h) <- NULL
  G <- vector("list", nrow(cross.mat))
  for(i in 1:nrow(cross.mat)){
    R <- array(NA, dim = c(nrow(h), mean(ploidy[cross.mat[i,]]), n.ind[i]),
               dimnames = list(h$mrks,NULL, paste0("Ind_",paste0(cross.mat[i,],collapse="x"), "_", 1:n.ind[i])))
    pid <- str_remove_all(cross.mat[i,], "P")
    Pa <- str_detect(str_split_fixed(h$mrks, "_", 2)[,1], pid[1])
    Pb <- str_detect(str_split_fixed(h$mrks, "_", 2)[,1], pid[2])
    Pab <- which(Pa & Pb)
    Ha <- str_detect(str_split_fixed(colnames(h), "_", 2)[,1], pid[1])
    Hb <- str_detect(str_split_fixed(colnames(h), "_", 2)[,1], pid[2])
    A <- simulate_cross(ploidy1 = ploidy[cross.mat[i,1]],
                        ploidy2 = ploidy[cross.mat[i,2]],
                        n.ind = n.ind[i],
                        n.mrk = length(Pab),
                        h1 = h[Pab,Ha],
                        h2 = h[Pab,Hb],
                        cm.map  = h[Pab,"map"])$offspring
    id <- apply(A, 1, function(x) length(unique(x)) > 1)
    dimnames(A) <- list(h[Pab,"mrks"], NULL, paste0("Ind_",paste0(cross.mat[i,],collapse="x"), "_", 1:n.ind[i]))
    for(j in 1:n.ind[i])
       R[rownames(A[id,,j]),,j] <- A[id,,j]
    G[[i]] <- R
  }
  names(G) <- apply(cross.mat, 1, paste0, collapse = "x")
  map <- h$map
  names(map) <- h$mrks
  ind.names <- unlist(lapply(G, function(x) dimnames(x)[[3]]))
  ped <- data.frame(Parent1 = str_split_fixed(ind.names, pattern = "_|x", 4)[,2],
                    Parent2 = str_split_fixed(ind.names, pattern = "_|x", 4)[,3],
                    row.names = ind.names)
  idp <- cbind(c(1, (cumsum(ploidy)+1)[-length(ploidy)]), cumsum(ploidy))
  H <- vector("list", length(ploidy))
  h1 <- h[,-c(1:2)]
  for(i in 1:length(H))
    H[[i]] <- data.frame(h1[,idp[i,1]:idp[i,2]], row.names = h$mrks)
  names(H) <- paste0("P", 1:length(H))
  structure(list(dat = G, ploidy = ploidy,
                 phases = H, pedigree = ped,
                 map = map),
            class = "mappoly2.data")
}

#' @rdname simulate_multiple_crosses
#' @importFrom reshape2 melt
#' @importFrom dplyr group_by summarise filter arrange n
#' @importFrom stringr %>%
#' @export
print.mappoly2.data <- function(x){
  y <- table(x$pedigree$Parent1, x$pedigree$Parent2)
  fds <- length(unique(unlist(dimnames(y))))
  n.mrk <- dim(x$dat[[1]])[[1]]
  w <- melt(x$dat)
  w <- w[!is.na(w$value),]
  w <- w %>% group_by(Var1) %>%
    summarise(cod = unique(value), .groups = 'drop') %>%
    group_by(Var1) %>%
    summarise(n = n(), .groups = 'drop') %>%
    arrange(desc(n))
  cat("This is an object of class mappoly2.data'\n")
  cat("    Number of founders:                     ", fds, "\n")
  cat("    Ploidy of founders:                     ", x$ploidy, "\n")
  cat("    No. individuals:                        ", sum(y), "\n")
  cat("    No. markers:                            ", n.mrk, "\n\n")
  cat("Number of individuals per crosses:")
  print(y)
  cat("\nNumber of markers per allelic level:")
  print(table(w$n))
}


drsimonj_colors <- c(
  `red`        = "#d11141",
  `green`      = "#00b159",
  `blue`       = "#00aedb",
  `orange`     = "#f37735",
  `yellow`     = "#ffc425",
  `light grey` = "#cccccc",
  `dark grey`  = "#8c8c8c")

my_pal <- c(darkslategray = "#2f4f4f",
            forestgreenmm = "#228b22",
            maroon2 = "#7f0000",
            indigo = "#4b0082",
            darkorange = "#ff8c00",
            yellow = "#ffff00",
            lime = "#00ff00",
            aqua = "#00ffff",
            blue = "#0000ff",
            cornflower = "#6495ed",
            moccasin = "#ffe4b5",
            hotpink = "#ff69b4")

#' Convert biallelic, single ploidy mappolymp data into mappoly legacy data
#'
#' @param void internal function to be documented
#' @examples
#' p <- c(4, 4, 4, 4)
#' names(p) <- c("P1", "P2", "P3", "P4")
#' cm <- matrix(c("P1","P2","P1","P3","P2","P3","P1","P1","P2","P4","P3","P4"),
#'              ncol = 2, byrow = T)
#' D <- simulate_multiple_crosses(ploidy = p, #four parents
#'                                cross.mat = cm,
#'                                n.ind = c(20, 20, 20, 20, 20, 20), # per cross
#'                                n.mrk= c(10,10,10,10), # per parent
#'                                alleles = list(c(0:1), c(0:1), c(0:1), c(0:1)),
#'                                map.length = 100)
#' D
#' bipar.dat <- mappolymp_to_mappoly(D)
#' require(mappoly)
#' plot(bipar.dat[[1]])
#'
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
mappolymp_to_mappoly <- function(dat){
  pop.table <- table(dat$pedigree[,1:2])
  pop.id <- as.logical(pop.table)
  pop.names <- apply(Reduce(expand.grid, dimnames(pop.table)), 1, paste0, collapse = "x")
  bipar.pops <- vector("list", sum(pop.id))
  names(bipar.pops) <- pop.names[pop.id]
  geno.dose <- sapply(dat$offspring, function(x) apply(x, 1, sum))
  for(i in 1:length(bipar.pops)){
    P1 <- strsplit(names(bipar.pops), "x")[[i]][1]
    P2 <- strsplit(names(bipar.pops), "x")[[i]][2]
    pl.ind <- apply(dat$pedigree[dat$pedigree[,1] == P1 & dat$pedigree[,2] == P2,3:4], 1, sum)/2
    pl <- unique(pl.ind)
    if(length(pl) > 1) stop("MAPpoly does no support mixed ploidies")
    geno.dose[is.na(geno.dose)] <- pl  + 1
    dosage.p1 <- apply(dat$phases[[P1]], 1, sum)
    dosage.p2 <- apply(dat$phases[[P2]], 1, sum)
    id.mrk.names <- names(which(!is.na(dosage.p1 + dosage.p2)))
    chrom <- rep(1, length(id.mrk.names))
    pos <- dat$map[dat$map$mrks%in%id.mrk.names,2]
    names(chrom) <- names(pos) <- id.mrk.names
    id.mrk.names <- id.mrk.names[!(dosage.p1[id.mrk.names] == 0 & dosage.p2[id.mrk.names] == 0 |
                                   dosage.p1[id.mrk.names] == pl & dosage.p2[id.mrk.names] == pl |
                                     dosage.p1[id.mrk.names] == pl & dosage.p2[id.mrk.names] == 0 |
                                     dosage.p1[id.mrk.names] == 0 & dosage.p2[id.mrk.names] == pl)]
    res <- structure(list(ploidy = pl,
                          n.ind = length(pl.ind),
                          n.mrk = length(id.mrk.names),
                          ind.names = names(pl.ind),
                          mrk.names = id.mrk.names,
                          dosage.p1 = dosage.p1[id.mrk.names],
                          dosage.p2 = dosage.p2[id.mrk.names],
                          chrom = chrom,
                          genome.pos = pos,
                          prob.thres = NULL,
                          geno.dose = geno.dose[id.mrk.names,names(pl.ind)],
                          nphen = 0,
                          phen = NULL,
                          chisq.pval = NULL),
                     class = "mappoly.data")
    res <- filter_non_conforming_classes(res)
    ##Computing chi-square p.values
    Ds <- array(NA, dim = c(pl+1, pl+1, pl+1))
    for(k in 0:pl)
      for(j in 0:pl)
        Ds[k+1,j+1,] <- segreg_poly(pl = pl, dP = k, dQ = j)
    Dpop <- cbind(res$dosage.p1, res$dosage.p2)
    M <- t(apply(Dpop, 1, function(x) Ds[x[1]+1, x[2]+1,]))
    dimnames(M) <- list(res$mrk.names, c(0:pl))
    M <- cbind(M, res$geno.dose)
    res$chisq.pval <- apply(M, 1, mrk_chisq_test, ploidy = pl)
    bipar.pops[[i]] <- res
  }
  return(bipar.pops)
}

#' mp2 to mp
#'
#' @param void internal function to be documented
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export mp2_csv_to_mappoly
mp2_csv_to_mappoly <- function(dat){
  ## Removing markers with missing data points for parents
  dat <- dat[which(!is.na(dat[,2,drop = TRUE]) & !is.na(dat[,3,drop = TRUE])),]
  ## get number of individuals -------------
  n.ind <- ncol(dat) - 5
  ## get number of markers -----------------
  n.mrk <- nrow(dat)
  ## get marker names ----------------------
  mrk.names <- dat[,1,drop = TRUE]
  ## get individual names ------------------
  ind.names <- colnames(dat)[-c(1:5)]
  ## get dosage in parent P ----------------
  dosage.p1 <- as.integer(dat[,2,drop = TRUE])
  ## get dosage in parent Q ----------------
  dosage.p2 <- as.integer(dat[,3,drop = TRUE])
  ## monomorphic markers
  dp <- abs(abs(dosage.p1-(ploidy/2))-(ploidy/2))
  dq <- abs(abs(dosage.p2-(ploidy/2))-(ploidy/2))
  #id <- dp + dq != 0
  ## get chromosome info ---------------------
  chrom <- as.character(dat[,4,drop = TRUE])
  ## get sequence position info ------------
  sequencepos <- as.numeric(dat[,5,drop = TRUE])
  names(sequencepos) <- names(chrom) <- names(dosage.p2) <- names(dosage.p1) <-  mrk.names
  nphen <- 0
  phen <- NULL
  ## get genotypic info --------------------
  geno.dose <- as.matrix(dat[,-c(1:5), drop = FALSE])
  dimnames(geno.dose) <- list(mrk.names, ind.names)
  geno.dose[is.na(geno.dose)] <- ploidy + 1
  ## returning the 'mappoly.data' object
  res <- structure(list(ploidy = ploidy,
                        n.ind = n.ind,
                        n.mrk = nrow(dat),
                        ind.names = ind.names,
                        mrk.names = mrk.names,
                        dosage.p1 = dosage.p1,
                        dosage.p2 = dosage.p2,
                        chrom = chrom,
                        genome.pos = sequencepos,
                        seq.ref = NULL,
                        seq.alt = NULL,
                        all.mrk.depth = NULL,
                        prob.thres = NULL,
                        geno.dose = geno.dose,
                        nphen = nphen,
                        phen = phen,
                        kept = NULL,
                        elim.correspondence = NULL),
                   class = "mappoly.data")
  return(res)
}

#' @export
get_full_sibs <- function(pedigree){
  pedigree <- pedigree[!apply(pedigree[,2:3], 1, function(x) any(is.na(x))), ]
  pedigree <- cbind(pedigree, Pop = apply(pedigree[, 2:3], 1, paste0, collapse = "-"))
  lapply(split(pedigree, pedigree$Pop),
         function(x){
           list(P1 = x[1,2],
                P2 = x[1,3],
                Ind  = x[,1])
         })
}

#' @export
create_map_multi_fam <- function(input.map, step = 0)
{
  map <- c(0, cumsum(imf_h(input.map$rf)))
  names(map) <- names(input.map$states)
  if(round(step, 1)  ==  0)
    return(map)
  minloc <- min(map)
  map <- map-minloc
  a <- seq(floor(min(map)), max(map), by = step)
  a <- a[is.na(match(a,map))]
  names(a) <- paste("loc",a,sep = "_")
  return(sort(c(a,map))+minloc)
}



#' Chi-square test
#'
#' @param void internal function to be documented
#' @keywords internal
mrk_chisq_test <- function(x, ploidy){
  y <- x[-c(1:(ploidy+1))]
  y[y == ploidy+1] <- NA
  y <- table(y, useNA = "always")
  names(y) <- c(names(y)[-length(y)], "NA")
  seg.exp <- x[0:(ploidy+1)]
  seg.exp <- seg.exp[seg.exp != 0]
  seg.obs <- seg.exp
  seg.obs[names(y)[-length(y)]] <- y[-length(y)]
  pval <- suppressWarnings(stats::chisq.test(x = seg.obs, p = seg.exp[names(seg.obs)])$p.value)
  pval
}


#' @export
#' @importFrom mappoly extract_map
plot_map_list_consensus <- function(map.list=NULL, consensus.map = NULL, col = "lightgray"){
  if(is.null(consensus.map))
  {
    if(length(col) == 1) col = rep(col, length(map.list))
    m <- vector("list", length(map.list))
    for(i in 1:length(map.list))
      m[[i]] <- extract_map(map.list[[i]])
    names(m) <- c(names(map.list))
    plot(0, ylim = c(0,length(m)+1), xlim = c(0,max(sapply(m, max))),
         type = "n",
         axes = FALSE,
         ylab = "", xlab = "Distance (cM)")
    axis(1)
    for(i in 1:length(m))
      mappoly:::plot_one_map(m[[i]], i, horiz = T, col = col[i])
    axis(2, at = 1:length(m), labels = names(m),
         lwd = 0, las = 2)
  }
  else if(is.null(map.list)){
    m <- vector("list", 1)
    m[[1]] <- cumsum(imf_h(c(0, consensus.map$rf)))
    names(m) <- c("cons")
    plot(0, ylim = c(0,length(m)+1), xlim = c(0,max(sapply(m, max))),
         type = "n",
         axes = FALSE,
         ylab = "", xlab = "Distance (cM)")
    axis(1)
    for(i in 1:length(m))
      mappoly:::plot_one_map(m[[i]], i, horiz = T, col = col[i])
    axis(2, at = 1:length(m), labels = names(m),
         lwd = 0, las = 2)
  }
  else{
    if(length(col) == 1) col = rep(col, length(map.list) + 1)
    m <- vector("list", length(map.list) + 1)
    m[[1]] <- cumsum(imf_h(c(0, consensus.map$rf)))
    for(i in 1:length(map.list))
      m[[i+1]] <- extract_map(map.list[[i]])
    names(m) <- c("cons", names(map.list))
    plot(0, ylim = c(0,length(m)+1), xlim = c(0,max(sapply(m, max))),
         type = "n",
         axes = FALSE,
         ylab = "", xlab = "Distance (cM)")
    axis(1)
    for(i in 1:length(m))
      mappoly:::plot_one_map(m[[i]], i, horiz = T, col = col[i])
    axis(2, at = 1:length(m), labels = names(m),
         lwd = 0, las = 2)
  }
}

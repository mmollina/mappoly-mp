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

#' Convert biallelic, single ploidy mappoly2 data into mappoly legacy data
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
#' bipar.dat <- mappoly2_to_mappoly(D)
#' require(mappoly)
#' plot(bipar.dat[[1]])
#'
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
mappoly2_to_mappoly <- function(dat){
  bipar.pops <- vector("list", length(dat$dat))
  names(bipar.pops) <- names(dat$dat)
  for(i in 1:length(bipar.pops)){
    P1 <- strsplit(names(dat$dat[i]), "x")[[1]][1]
    P2 <- strsplit(names(dat$dat[i]), "x")[[1]][2]
    geno.dose <- apply(dat$dat[[i]], c(1,3), sum)
    geno.dose[is.na(geno.dose)] <- dat$ploidy + 1
    dosage.p1 <- apply(dat$phases[[P1]], 1, sum)
    dosage.p2 <- apply(dat$phases[[P2]], 1, sum)
    id.mrk.names <- names(which(!is.na(dosage.p1 + dosage.p2)))
    bipar.pops[[i]] <- structure(list(ploidy = dat$ploidy[P1],
                                      n.ind = dim(dat$dat[[i]])[3],
                                      n.mrk = length(id.mrk.names),
                                      ind.names = dimnames(dat$dat[[i]])[[3]],
                                      mrk.names = id.mrk.names,
                                      dosage.p1 = dosage.p1[id.mrk.names],
                                      dosage.p2 = dosage.p2[id.mrk.names],
                                      chrom = rep(1, length(id.mrk.names)),
                                      genome.pos = dat$map[id.mrk.names],
                                      prob.thres = NULL,
                                      geno.dose = geno.dose[id.mrk.names,],
                                      nphen = 0,
                                      phen = NULL,
                                      chisq.pval = NULL),
                                 class = "mappoly.data")
  }
  return(bipar.pops)
}





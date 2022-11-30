#' Prepare maps to merge
#'
#' @param map.list a list of maps to merge. All objects must be of class
#'      \code{mappoly.map}
#' @param data.list a list of objects of class \code{mappoly.data} corresponding
#'      to the maps in argument \code{map.list}
#' @param parents.mat a matrix where each row contains the names of the parents
#'      corresponding to the maps in argument \code{map.list}
#'@examples
#'    \donttest{
#'    map_ch1_BExMG <- readRDS("~/repos/current_work/rose/fullsib_maps/BExMG/map_err_ch_1.rds")
#'    map_ch1_SWxBE <- readRDS("~/repos/current_work/rose/fullsib_maps/SWxBE/map_err_ch_1.rds")
#'
#'
#'
#'    map.list <- list(BExMG = map_ch1_BExMG,
#'                     SWxBE = map_ch1_SWxBE)
#'    dat_BExMG <- readRDS("~/repos/current_work/rose/data/dat_BExMG.rds")
#'    dat_SWxBE <- readRDS("~/repos/current_work/rose/data/dat_SWxBE.rds")
#'    data.list <- list(BExMG = dat_BExMG,
#'                      SWxBE = dat_SWxBE)
#'    dat <- dat_BExMG
#'    map_ch1_BExMG <- get_submap(map_ch1_BExMG, mrk.pos = 1:100, reestimate.rf = FALSE)
#'    dat <- dat_SWxBE
#'    map_ch1_SWxBE <- get_submap(map_ch1_SWxBE, mrk.pos = 1:100, reestimate.rf = FALSE)
#'    map.list <- list(BExMG = map_ch1_BExMG,
#'                     SWxBE = map_ch1_SWxBE)
#'    parents.mat <- matrix(c("BE", "MG", "SW", "BE"),
#'                          2, 2, byrow = TRUE,
#'                          dimnames = list(c("pop1", "pop2"), c("P1", "P2")))
#'   }
#'
#' @importFrom plotly ggplotly
#' @export
prepare_maps_to_merge <- function(map.list,
                                  data.list,
                                  parents.mat){
  pl <- data.list[[1]]$ploidy
  ## Checking input arguments
  if (any(!sapply(map.list, inherits, "mappoly.map")))
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if (any(!sapply(data.list, inherits, "mappoly.data")))
    stop("All elemnts in 'data.list' should be of class 'mappoly.data'")
  ## Gathering parent's phases
  w <- table(as.vector(parents.mat))
  phases <- vector("list", length(w))
  names(phases) <- names(w)
  for(i in 1:length(phases)){
      par.ord <- which(parents.mat == names(phases)[i], arr.ind = T)
      hom.res <- match_homologs(map.list, par.ord)
      phases[[i]] <- apply(hom.res$ph, 1, paste, collapse = "")
  }
  ## Gathering pedigree
  pedigree <- NULL
  for(i in 1:nrow(parents.mat)){
    pedigree<- rbind(pedigree,
                     data.frame(Off = data.list[[i]]$ind.names,
                                Par1 = parents.mat[i,1],
                                Par2 = parents.mat[i,2],
                                pl1 = data.list[[i]]$ploidy,
                                pl2 = data.list[[i]]$ploidy,
                                pop = i,
                                row.names = data.list[[i]]$ind.names))
  }
  pedigree <- pedigree[,-1]
  ## Gathering offspring genotypes
  R <- generate_biallelic_indices(pl)
  emit <- states <- vector("list", length(phases[[1]]))
  names(emit) <- names(states) <- names(phases[[1]])
  for(j in names(phases[[1]])){
    Etemp <- Ltemp <- vector("list", nrow(pedigree))
    names(Etemp) <- names(Ltemp) <- rownames(pedigree)
    for(i in rownames(pedigree)){
      ##FIXME Split in full-sibs
      id<- paste(phases[[pedigree[i,"Par1"]]][j],
                 phases[[pedigree[i,"Par2"]]][j],
                 sep = "-")
      if(str_detect(id, "NA"))
        Ltemp[[i]] <- R[[1]][[1]]
      else{
        x <- data.list[[pedigree[i, "pop"]]]$geno.dose[j,i]
        if(x == pl + 1)
          Ltemp[[i]] <- R[[1]][[1]]
        else
          Ltemp[[i]] <- R[[id]][[as.character(x)]]
      }
      Etemp[[i]] <- matrix(rep(1, nrow(Ltemp[[i]])), ncol = 1)
    }
    states[[j]] <- Ltemp
    emit[[j]] <- Etemp
  }
  list(n.mrk = length(phases[[1]]),
       n.ind = nrow(pedigree),
       states = states,
       emit = emit,
       ploidy = pedigree[,c("pl1", "pl2")],
       genome.pos = hom.res)
}

match_homologs <- function(map.list, par.ord){
  #### Check order of parents
  ## Shared markers
  idn <- Reduce(intersect, lapply(map.list, function(x) x$info$mrk.names))
  ## Markers from both maps and yheir positions
  id.all <- unlist(lapply(map.list, function(x) x$info$mrk.names))
  pos.all <- unlist(lapply(map.list, function(x) x$info$genome.pos))
  ## Removing duplicate markers
  w <- unique(data.frame(id.all=id.all, pos.all=pos.all, row.names = NULL))
  ## Ordering according genome
  pos <- w[order(w$pos.all),]
  ## Gathering phases
  ph.list <- ph.mat <- NULL
  for(i in 1:nrow(par.ord)){
    pl <- map.list[[par.ord[i,1]]]$info$ploidy
    ph <- map.list[[par.ord[i,1]]]$maps[[1]]$seq.ph[[par.ord[i,2]]]
    ph <- ph_list_to_matrix(ph, pl)
    dimnames(ph) <- list(map.list[[par.ord[i,1]]]$info$mrk.names, paste0(letters[1:pl], i))
    ph.list[[i]] <- ph
    ph.mat <- rbind(ph.mat, t(ph[idn,]))
  }
  if(nrow(par.ord) > 1){
    dd <- dist(ph.mat, method = "binary")
    hc <- hclust(dd, method = "ward.D2")
    d <- as.dendrogram(hc)
    my_pal <- c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07")
    d <- d %>%
      dendextend::color_branches(k=4, col = my_pal) %>%
      dendextend::color_labels(k = 4, col = my_pal)
    plot(d)
    dendextend::rect.dendrogram(d, k = 4, lwd = 3, border = my_pal)
    homologs  <- cutree(hc, k = 4)
    y <- split(names(homologs), as.factor(homologs))
    names(y) <- paste0("h", 1:pl)
    x <- y %>%
      melt %>%
      mutate(homolog = substr(value, 1,1)) %>%
      mutate(pop = substr(value,2,10)) %>%
      arrange(pop, L1) %>%
      acast(pop ~ L1, value.var = "value")
    ## Re-organizing homologs
    for(i in 1:length(ph.list))
      ph.list[[i]] <- ph.list[[i]][,x[i,]]
    ##FIXME: Instead of selecting the first one,
    ## use some probability or threshold
    id.ph <- logical(length(idn))
    for(i in 1:length(idn))
      id.ph[i] <- as.logical(length(unique(unlist(lapply(ph.list, function(x) paste0(x[1, ], collapse = ""))))))
    ph.out <- ph.list[[1]][idn[id.ph],]
    for(i in 1:length(ph.list)){
      idtemp <- setdiff(rownames(ph.list[[i]]), rownames(ph.out))
      ph.out <- rbind(ph.out, ph.list[[i]][idtemp,])
    }
    if(length(idn[!id.ph]) > 0)
      ph.out[idn[!id.ph],][] <- NA
  } else {
    ph.out <- ph.list[[1]]
    for(i in 1:length(map.list)){
      idtemp <- setdiff(map.list[[i]]$info$mrk.names, rownames(ph.out))
      ph.out <- rbind(ph.out, matrix(NA, length(idtemp), pl, dimnames = list(idtemp, NULL)))
      x <- d <- NULL
    }
  }
  ph.out <- ph.out[pos$id.all,]
  colnames(pos) <- c("mrk", "geno.pos")
  list(ph = ph.out,
       equivalence = x,
       dend = d,
       shared.mrks = idn,
       genome.pos = pos)
}

generate_biallelic_indices <- function(pl){
  I <- matrix(0, 1, pl)
  for(d in 1:pl){
    x <- combn(1:pl, d)
    y <- matrix(0, ncol(x), pl)
    for(i in 1:ncol(x))
      y[i,x[,i]][] <- 1
    I <- rbind(I, y)
  }
  I <- rbind(I, matrix(1, 1, pl))
  Z <- vector("list", nrow(I)^2)
  n <- character(nrow(I)^2)
  cte <- 1
  a <- apply(I, 1, paste, collapse="")
  A <- as.matrix(expand.grid(1:choose(pl, pl/2), 1:choose(pl,pl/2))[,2:1] - 1)
  for(i in 1:nrow(I)){
    for(j in 1:nrow(I)){
      n[cte] <- paste(a[i], a[j], sep = "-")
      v1 <- apply(combn(I[i,], pl/2), 2, sum)
      v2 <- apply(combn(I[j,], pl/2), 2, sum)
      v <- kronecker(v1, v2, "+")
      u <- order(v)
      u <- split(u, as.factor(sort(v)))
      for(k in 1:length(u))
        u[[k]] <- A[u[[k]], , drop = FALSE]
      Z[[cte]] <- u
      cte <- cte + 1
    }
  }
  names(Z) <- n
  Z
}

compare_single_maps <- function(map.list){
  mrk.id <- unlist(lapply(map.list, function(x) x$info$mrk.names))
  geno_pos <- unlist(lapply(map.list, function(x) x$info$genome.pos))
  L.temp <- unique(data.frame(mrk.id = mrk.id,
                         geno_pos = geno_pos,
                         sapply(map.list, function(x, i) extract_map(x)[i], i = mrk.id),
                         row.names = NULL))
  rg0 <- range(L.temp[,-c(1:2)], na.rm = T)
  L.temp$u.s<-apply(L.temp[,-c(1:2)], 1, function(x) {
    if(is.na(x[1]))
      return(0)
    else if(is.na(x[2]))
      return(1)
    else
      return(2)
  })

  L.temp$colorVal <- as.factor(L.temp$u.s)
  levels(L.temp$colorVal) <- 1:length(levels(L.temp$colorVal))

  dm<- vector("list", ncol(L.temp)-3)
  for(i in 2:(length(dm)+1)){
    if(colnames(L.temp)[i] == "geno_pos")
      rg <- range(L.temp[,"geno_pos"], na.rm = T)
    else rg <- rg0
    dm[[i]] <-     list(range = rg,
                        label = colnames(L.temp)[i],
                        values = as.formula(paste0("~",colnames(L.temp)[i])))
  }
  figA <- L.temp %>% plot_ly(type = 'parcoords',
                             line = list(color = ~colorVal,
                                         lwd = .2,
                                         colorscale = list(c(0, "#B2DF8A"),
                                                           c(0.5, "#FB9A99"),
                                                           c(1, "#1F78B4"))),
                             dimensions = dm
  )
  figA
}

dose_2_vec<- function(d, pl){
  if(is.na(d))
    return(rep(NA, pl))
  if(d == 0)
    return(rep(0, pl))
  if(d > pl | d < 0)
    return(rep(NA, pl))
  x<- rep(0, pl)
  x[1:d] <- 1
  return(x)
}

# prepare_maps_to_merge2 <- function(map.list,
#                                    data.list,
#                                    parents.mat){
#   pl <- data.list[[1]]$ploidy
#   ## Checking input arguments
#   if (any(!sapply(map.list, inherits, "mappoly.map")))
#     stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
#   if (any(!sapply(data.list, inherits, "mappoly.data")))
#     stop("All elemnts in 'data.list' should be of class 'mappoly.data'")
#   ## Gathering parent's phases
#   w <- table(as.vector(parents.mat))
#   phases <- vector("list", length(w))
#   names(phases) <- names(w)
#   for(i in 1:length(phases)){
#     par.ord <- which(parents.mat == names(phases)[i], arr.ind = T)
#     hom.res <- match_homologs(map.list, par.ord)
#     phases[[i]] <- hom.res$ph
#     #phases[[i]] <- apply(hom.res$ph, 1, paste, collapse = "")
#   }
#   #lapply(phases, head)
#   ## Gathering pedigree
#   pedigree <- NULL
#   for(i in 1:nrow(parents.mat)){
#     pedigree<- rbind(pedigree,
#                      data.frame(Off = data.list[[i]]$ind.names,
#                                 Par1 = parents.mat[i,1],
#                                 Par2 = parents.mat[i,2],
#                                 pl1 = data.list[[i]]$ploidy,
#                                 pl2 = data.list[[i]]$ploidy,
#                                 row.names = data.list[[i]]$ind.names))
#   }
#   pedigree <- pedigree[,-1]
#   ## Gathering offspring genotypes
#   L <- vector("list", length(map.list))
#   for(i in 1:length(L)){
#     L[[i]] <- as.data.frame(data.list[[i]]$geno.dose[map.list[[i]]$info$mrk.names, ]) %>%
#       rownames_to_column("mrk.id")
#   }
#   D <-purrr::reduce(L, dplyr::full_join, by = "mrk.id") %>%
#     column_to_rownames("mrk.id")
#   D <- D[rownames(hom.res$ph), ]
#   #image(as.matrix(D))
#   offspring <- vector("list", ncol(D))
#   names(offspring) <- colnames(D)
#   for(i in 1:ncol(D))
#     offspring[[i]] <- t(apply(D[,i,drop = F], 1, dose_2_vec, pl = pl))
#   w <- list(offspring = offspring,
#             phases = phases,
#             pedigree = pedigree,
#             map = NA,
#             joint.info = NA)
#   states <- states_to_visit(w)
#   states
# }
#

#' Prepare maps to merge
#'
#' @param map.list a list of maps to merge. All objects must be of class
#'      \code{mappoly.map}
#' @param data.list a list of objects of class \code{mappoly.data} corresponding
#'      to the maps in argument \code{map.list}
#' @param parents.mat a matrix where each row contains the names of the parents
#'      corresponding to the maps in argument \code{map.list}
#' @param err global genotyping error
#'
#'@examples
#'    \donttest{
#'    map_ch1_BExMG <- readRDS("~/repos/current_work/rose/fullsib_maps/BExMG/map_err_ch_1.rds")
#'    map_ch1_SWxBE <- readRDS("~/repos/current_work/rose/fullsib_maps/SWxBE/map_err_ch_1.rds")
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
#' @export
prepare_maps_to_merge <- function(map.list,
                                  data.list,
                                  parents.mat,
                                  err = 0){
  ## Checking input arguments
  if (any(!sapply(map.list, inherits, "mappoly.map")))
    stop("All elemnts in 'map.list' should be of class 'mappoly.map'")
  if (any(!sapply(data.list, inherits, "mappoly.data")))
    stop("All elemnts in 'data.list' should be of class 'mappoly.data'")

  ## Gathering parent's phases
  w <- table(as.vector(parents.mat))
  hom.res <- phases <- vector("list", length(w))
  names(phases) <- names(w)

  ## Gathering ploidy level
  pl <- data.list[[1]]$ploidy

  ## For each unique parent
  for(i in names(phases)){
      par.ord <- which(parents.mat == i, arr.ind = T)
      hom.res[[i]] <- match_homologs(map.list, par.ord, pl)
      phases[[i]] <- apply(hom.res[[i]]$ph, 1, paste, collapse = "")
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
  E <- R[[2]]
  R <- R[[1]]
  emit <- states <- vector("list", length(phases[[1]]))
  names(emit) <- names(states) <- names(phases[[1]])

  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = length(phases[[1]]), # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  cte <- 0
  for(j in names(phases[[1]])){
    Etemp <- Ltemp <- vector("list", nrow(pedigree))
    names(Etemp) <- names(Ltemp) <- rownames(pedigree)
    for(i in rownames(pedigree)){
      if(err <= 0.001) {
        if(!j%in%data.list[[pedigree[i, "pop"]]]$mrk.names){
          Ltemp[[i]] <- R[[1]][[1]]
          Etemp[[i]] <- matrix(rep(1, nrow(Ltemp[[i]])), ncol = 1)
          next()
        }
        if(!j%in%map.list[[pedigree[i, "pop"]]]$info$mrk.names){
          Ltemp[[i]] <- R[[1]][[1]]
          Etemp[[i]] <- matrix(rep(1, nrow(Ltemp[[i]])), ncol = 1)
          next()
        }
        ##FIXME Split in full-sibs
        id<- paste(phases[[pedigree[i,"Par1"]]][j],
                   phases[[pedigree[i,"Par2"]]][j],
                   sep = "-")
        if(str_detect(id, "NA"))
          Ltemp[[i]] <- R[[1]][[1]]
        else{
          x <- data.list[[pedigree[i, "pop"]]]$geno.dose[j,i]
          if(!x%in%names(R[[id]])){
            Ltemp[[i]] <- R[[1]][[1]]
            Etemp[[i]] <- matrix(rep(1, nrow(Ltemp[[i]])), ncol = 1)
            next()
          }
          if(x > pl)
            Ltemp[[i]] <- R[[1]][[1]]
          else
            Ltemp[[i]] <- R[[id]][[as.character(x)]]
        }
        if(is.null(Ltemp[[i]]))
          Ltemp[[i]] <- R[[1]][[1]]
        Etemp[[i]] <- matrix(rep(1, nrow(Ltemp[[i]])), ncol = 1)
      }
      else if (err > 0.001) {
        Ltemp[[i]] <- R[[1]][[1]]
        if(!j%in%data.list[[pedigree[i, "pop"]]]$mrk.names){
          Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                   nrow(Ltemp[[i]])),
                               ncol = 1)
          next()
        }
        if(!j%in%map.list[[pedigree[i, "pop"]]]$info$mrk.names){
          Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                   nrow(Ltemp[[i]])),
                               ncol = 1)
          next()
        }
        ##FIXME Split in full-sibs
        id<- paste(phases[[pedigree[i,"Par1"]]][j],
                   phases[[pedigree[i,"Par2"]]][j],
                   sep = "-")
        if(str_detect(id, "NA"))
          Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                   nrow(Ltemp[[i]])),
                               ncol = 1)
        else{
          x <- data.list[[pedigree[i, "pop"]]]$geno.dose[j,i]
          if(!x%in%names(R[[id]])){
            Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                     nrow(Ltemp[[i]])),
                                 ncol = 1)
            next()
          }
          if(x == pl + 1)
            Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                     nrow(Ltemp[[i]])),
                                 ncol = 1)
          else{
            v <- E[[id]][[as.character(x)]]
            A <- matrix(rep(0, nrow(Ltemp[[i]])),
                                 ncol = 1)
            A[v,1] <- rep((1-err)/length(v))
            A[A[,1] == 0,1] <- err/sum(A[,1] == 0)
            Etemp[[i]] <- A
          }
        }
        if(is.null(Ltemp[[i]])){
          Etemp[[i]] <- matrix(rep(1/nrow(Ltemp[[i]]),
                                   nrow(Ltemp[[i]])),
                               ncol = 1)
        }
      } else stop("should not get here!")
    }
    states[[j]] <- Ltemp
    emit[[j]] <- Etemp
    cte <- cte + 1
    setTxtProgressBar(pb,  cte)
  }
  close(pb)
  structure(list(n.mrk = length(phases[[1]]),
       n.ind = nrow(pedigree),
       states = states,
       emit = emit,
       ploidy = pedigree[,c("pl1", "pl2")],
       merged.phases = hom.res,
       pedigree = pedigree,
       err = err), class = "multi.mappoly.consensus.info")
}

#' @export
plot.multi.mappoly.consensus.info <- function(w){
  hc <- sapply(w$merged.phases, function(x) x$hc)
  hc <-  hc[!sapply(hc, function(x) all(is.na(x)))]
  u1 <- unique(w$pedigree[,c(1,3)])
  v1 <- u1[,2]
  names(v1) <- u1[,1]
  u2 <- unique(w$pedigree[,c(2,4)])
  v2 <- u2[,2]
  names(v2) <- u2[,1]
  v <- c(v1, v2)
  a<-sqrt(length(hc))
  par(mfrow = c(floor(a), ceiling(a)))
  for(i in 1:length(hc)){
    d <- as.dendrogram(hc[[i]])
    d <- d %>%
      dendextend::color_branches(k = v[names(hc)[1]], col = mp_pallet2(v[names(hc)[1]])) %>%
      dendextend::color_labels(k = v[names(hc)[1]], col = mp_pallet2(v[names(hc)[1]]))
    plot(d, main = names(hc)[i])
    dendextend::rect.dendrogram(d, k = v[names(hc)[1]], lwd = 3,
                                border = mp_pallet2(v[names(hc)[1]]))
  }
}

#' @importFrom dendextend color_branches color_labels rect.dendrogram
#' @importFrom reshape2 acast
#' @export
match_homologs <- function(map.list, par.ord, pl){
  ## Number of full-sibs with the analyzed parent
  n <- nrow(par.ord)
  ## Shared markers
  idn <- Reduce(intersect, lapply(map.list, function(x) x$info$mrk.names))

  ## Markers from all maps and their positions
  id.all <- unlist(lapply(map.list, function(x) x$info$mrk.names))
  pos.all <- unlist(lapply(map.list, function(x) x$info$genome.pos))

  ## Removing duplicate markers
  w <- unique(data.frame(id.all=id.all, pos.all=pos.all, row.names = NULL))
  ## Ordering according genome
  pos <- w[order(w$pos.all),]
  ## Gathering phases
  ph.list <- vector("list", n)
  ph.mat <- NULL
  for(i in 1:n){
    ph <- map.list[[par.ord[i,1]]]$maps[[1]]$seq.ph[[par.ord[i,2]]]
    ph <- ph_list_to_matrix(ph, pl)
    dimnames(ph) <- list(map.list[[par.ord[i,1]]]$info$mrk.names, paste0(letters[1:pl], i))
    ph.list[[i]] <- ph
    ph.mat <- rbind(ph.mat, t(ph[idn,]))
  }
  if(n > 1){
    dd <- as.matrix(dist(ph.mat, method = "binary"))
    for(i in 1:n){
      id <- (((i-1)*pl)+1):(pl*i)
      dd[id,id][] <- 1
    }
    #image(dd)
    dd <- as.dist(dd)
    hc <- hclust(dd, method = "ward.D2")
    homologs  <- cutree(hc, k = pl)
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

    #id<-Reduce(intersect, sapply(ph.list, rownames))
    #image(ph.list[[1]][id,])
    #image(ph.list[[2]][id,])
    #image(ph.list[[3]][id,])

    ## Check if dose and phase is the same across shared markers
    S <- sapply(ph.list, function(x) apply(x[idn,], 1, paste0, collapse = ""))
    id.ph <- which(apply(S, 1, function(x) length(unique(x)) != 1))
    ph.out <- ph.list[[1]]
    ph.out[names(id.ph),][] <- NA
    for(i in 1:n){
      idtemp <- setdiff(rownames(ph.list[[i]]), rownames(ph.out))
      ph.out <- rbind(ph.out, ph.list[[i]][idtemp,])
    }
    remaining <- setdiff(pos$id.all, rownames(ph.out))
    if(length(remaining) > 0)
      ph.out <- rbind(ph.out,
                      matrix(NA, length(remaining), pl, dimnames = list(remaining, colnames(ph.out))))
  }
  else {
    hc <- NA
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
       hc = hc,
       shared.mrks = idn,
       genome.pos = pos)
}

#' @export
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
  W <- Z <- vector("list", nrow(I)^2)
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
      W[[cte]] <- u <- split(u, as.factor(sort(v)))
      for(k in 1:length(u)){
        u[[k]] <- A[u[[k]], , drop = FALSE]
      }
      Z[[cte]] <- u
      cte <- cte + 1
    }
  }
  names(W) <- names(Z) <- n
  list(Z, W)
}

#' @export
#' @importFrom mappoly extract_map
#' @importFrom plotly plot_ly
compare_single_maps <- function(map.list){
  mrk.id <- unlist(lapply(map.list, function(x) x$info$mrk.names))
  geno_pos <- unlist(lapply(map.list, function(x) x$info$genome.pos))
  L.temp <- unique(data.frame(mrk.id = mrk.id,
                         geno_pos = geno_pos,
                         sapply(map.list, function(x, i) extract_map(x)[i], i = mrk.id),
                         row.names = NULL))
  rg0 <- range(L.temp[,-c(1:2)], na.rm = T)
  L.temp$u.s <- apply(L.temp[,-c(1:2)], 1, function(x) sum(!is.na(x)))
  L.temp$u.s <- round(L.temp$u.s/max(L.temp$u.s),2)
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
                                         colorscale = 'YlGnBu',
                                         showscale = TRUE), #all
                             dimensions = dm) %>%
    layout(plot_bgcolor='rgb(211, 211, 211)') %>%
    layout(paper_bgcolor='rgb(211, 211, 211)')
  figA
}

#' @export
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




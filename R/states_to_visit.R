#' Generates the states the model should visit including global error rate
#'
#' @param void internal function to be documented
#' @examples
#' require(tidyverse)
#' p <- c(4, 2, 6, 4)
#' names(p) <- c("P1", "P2", "P3", "P4")
#' cm <- matrix(c("P1","P2","P1","P3","P2","P3","P1","P1","P2","P4","P3","P4"),
#'              ncol = 2, byrow = T)
#' D <- simulate_multiple_crosses(ploidy = p, #four parents
#'                                cross.mat = cm,
#'                                n.ind = c(20, 20, 20, 20, 20, 20), # per cross
#'                                n.mrk= c(10,10,10,10), # per parent
#'                                alleles = list(c(0:3), c(4:5), c(6:11), c(12:15)),
#'                                map.length = 100)
#' D
#' H <- states_to_visit(D, err = 0.05, is.log = TRUE)
#' names(H)
#' head(H$hmm.info)
#' mrk.sample <- names(D$map)
#' i1 <- H$hmm.info %>%
#'       filter(ind == "Ind_P3xP4_1", mrk %in% mrk.sample) %>%
#'       mutate(st.all = paste(st.p1, st.p2, sep = "-")) %>%
#'       mutate(x = 1/emit)
#' i1$mrk <- as_factor(i1$mrk)
#' p <- ggplot(i1, aes(x=st.all, y = x, fill = emit)) +
#'      scale_fill_distiller(palette = "RdYlBu") +
#'      geom_bar(stat="identity") + facet_grid(.~mrk)
#' p + coord_flip()
#'
#' H <- states_to_visit(D, err = 0, is.log = FALSE)
#' head(H$hmm.info)
#' mrk.sample <- names(D$map)
#' i1 <- H$hmm.info %>%
#'       filter(ind == "Ind_P1xP2_1", mrk %in% mrk.sample) %>%
#'       mutate(st.all = paste(st.p1, st.p2, sep = "-"))
#' i1$mrk <- as_factor(i1$mrk)
#' p <- ggplot(i1, aes(x=st.all, y = emit, fill = emit)) +
#'      scale_fill_distiller(palette = "YlOrBr",direction = 1) +
#'      geom_bar(stat="identity") + facet_grid(.~mrk)
#' p + coord_flip()
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @importFrom dplyr bind_rows
#' @importFrom forcats as_factor
#' @export
states_to_visit <- function(input.data, err = 0.0, is.log = TRUE){
  ## parsing input object
  ploidy <- input.data$ploidy
  phases <- input.data$phases
  pedigree <- input.data$pedigree
  mrks <- names(input.data$map)
  dat <- input.data$dat
  ngam <- choose(ploidy, ploidy/2)
  crosses <- str_split_fixed(names(input.data$dat),"x", 2)
  st <- 1:nrow(crosses)
  pls <- cbind(ploidy[crosses[,1]], ploidy[crosses[,2]])
  nst <- apply(pls, 1, paste, collapse = "x")
  ost <- order(nst)
  pls <- pls[ost, ,drop = FALSE]
  rownames(pls) <- names(st) <- nst[ost]
  Y <- vector("list", sum(apply(sapply(dat, dim)[-2,,drop = FALSE],2,prod)))
  cte <- 1
  ## For each marker
  for(i in mrks){
    geno.cur.mrk <- lapply(dat, function(x, i) x[i,,], i)
    #for(k in 1:length(geno.cur.mrk))
     #  geno.cur.mrk[[k]] <- apply(geno.cur.mrk[[k]], 2, sort, decreasing = TRUE)
    ## For all crosses
    for(k in 1:nrow(crosses)){
      cross.ind.names <- colnames(geno.cur.mrk[[k]])
      p1 <- crosses[k,1]
      p2 <- crosses[k,2]
      P1M <- phases[[p1]]
      P2M <- phases[[p2]]
      S <- data.frame(expand.grid(0:(ngam[p2]-1),
                                  0:(ngam[p1]-1))[,2:1],
                      as.numeric(1/(ngam[p1] * ngam[p2])))
      if(any(is.na(c(unlist(P1M[i,]), unlist(P2M[i,]))))){ ## If no crosses for these parents
        for(j in cross.ind.names){
          X <- data.frame(ind = j,
                          mrk = i,
                          pl1 = ploidy[p1],
                          pl2 = ploidy[p2],
                          pl.id = st[paste(ploidy[p1], ploidy[p2], sep = "x")],
                          st.p1 = S[,1],
                          st.p2 = S[,2],
                          st.all = 1:nrow(S)-1,
                          emit = S[,3],
                          row.names = NULL)
          if(is.log) X$emit <- log(X$emit)
          Y[[cte]] <- X
          cte <- cte + 1
        }
      } else {
        A1 <- combn(as.numeric(P1M[i,]), ploidy[p1]/2)
        A2 <- combn(as.numeric(P2M[i,]), ploidy[p2]/2)
        w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                       apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
        ### Indexed states
        a1 <- apply(str_split_fixed(w,
                                    "_",
                                    ploidy[p1]/2 + ploidy[p2]/2),
                    1, function(x)
                      paste0(sort(as.numeric(x)),
                             collapse = "_"))
        ### Observed genotype
        d <- apply(geno.cur.mrk[[k]],
                   2, function(x)
                     paste0(sort(as.numeric(x)), collapse = "_"))
        for(j in cross.ind.names){
          if(is.na(d[j]) | !d[j]%in%a1){
            X <- data.frame(ind = j,
                            mrk = i,
                            pl1 = ploidy[p1],
                            pl2 = ploidy[p2],
                            pl.id = st[paste(ploidy[p1], ploidy[p2], sep = "x")],
                            st.p1 = S[,1],
                            st.p2 = S[,2],
                            st.all = 1:nrow(S)-1,
                            emit = S[,3],
                            row.names = NULL)
            if(is.log) X$emit <- log(X$emit)
          } else {
            if(err < 1e-4) {
              X <- data.frame(ind = j,
                              mrk = i,
                              pl1 = ploidy[p1],
                              pl2 = ploidy[p2],
                              pl.id = st[paste(ploidy[p1], ploidy[p2], sep = "x")],
                              st.p1 = S[a1==d[j], 1],
                              st.p2 = S[a1==d[j], 2],
                              st.all = which(a1==d[j])-1,
                              emit = S[a1==d[j], 3],
                              row.names = NULL)
              X$emit <- 1/nrow(X)
              if(is.log) X$emit <- log(X$emit)
            } else{
              X <- data.frame(ind = j,
                              mrk = i,
                              pl1 = ploidy[p1],
                              pl2 = ploidy[p2],
                              pl.id = st[paste(ploidy[p1], ploidy[p2], sep = "x")],
                              st.p1 = S[,1],
                              st.p2 = S[,2],
                              st.all = 1:nrow(S)-1,
                              emit = S[,3],
                              row.names = NULL)
              X$emit[a1==d[j]] <- (1-err)/sum(a1==d[j])
              X$emit[a1!=d[j]] <- (err)/sum(a1!=d[j])
              if(is.log) X$emit <- log(X$emit)
            }
          }
          if(is.null(X)){
            X <- data.frame(ind = j,
                            mrk = i,
                            pl.1 = ploidy[p1],
                            pl.2 = ploidy[p2],
                            pl.id = st[paste(ploidy[p1], ploidy[p2], sep = "x")],
                            st.p1 = S[,1],
                            st.p2 = S[,2],
                            st.all = 1:nrow(S)-1,
                            emit = S[,3],
                            row.names = NULL)
            if(is.log) X$emit <- log(X$emit)
          }
          Y[[cte]] <- X
          cte <- cte + 1
        }
      }
    }
  }
  Y <- Y %>% bind_rows %>% arrange(ind)
  Y$ind <- as_factor(Y$ind)
  Y$mrk <- as_factor(Y$mrk)
  list(hmm.info = Y, err = err, is.log = is.log, ploidy.cross.id = pls)
}

#' Generates the states the model should visit for mappoly legacy
#'
#' @param void internal function to be documented
#' @examples
#' require(mappoly)
#' n.mrk <- 2
#' h.temp <- sim_homologous(ploidy = 4,
#'                          n.mrk = n.mrk,
#'                          max.d = 2,
#'                          max.ph = 0,
#'                          seed = 1)
#'
#' for(i in 1:n.mrk){
#'   h.temp$hom.allele.p[[i]] <- 1:2
#'   h.temp$hom.allele.q[[i]] <- 1
#'   h.temp$p[i] <- 2
#'   h.temp$q[i] <- 1
#' }
#' dat <- poly_cross_simulate(ploidy = 4,
#'                            rf.vec = .01,
#'                            n.mrk = n.mrk,
#'                            n.ind = 10,
#'                            h.temp,
#'                            seed = 8532,
#'                            draw = TRUE)
#' sim.map<-cumsum(c(0,rep(imf_h(.01), (n.mrk - 1))))
#' plot(dat)
#' s <- make_seq_mappoly(dat, "all")
#' tpt <- est_pairwise_rf(s)
#' tpt$pairwise$`1-2`
#' map <- est_rf_hmm_sequential(input.seq = s, twopt = tpt)
#' plot(map)
#' mp <- round(cumsum(mappoly::imf_h(c(0, map$maps[[1]]$seq.rf))),2)
#' ph <- list(ph_list_to_matrix(h.temp$hom.allele.p, dat$ploidy),
#'            ph_list_to_matrix(h.temp$hom.allele.q, dat$ploidy))
#' states.hmm <- states_to_visit_mp1(dat, ph, is.log = TRUE)
#' x1 <- est_map_R(states.hmm,tol = 1e-3, verbose = FALSE)
#' mp2 <- round(cumsum(mappoly::imf_h(c(0, x1[[2]]))), 2)
#' Y1 <- rbind(sim.map,mp,mp2)
#' rownames(Y1) <- c("simulation", "mappoly", "mappoly2")
#' Y1
#'
#' @author Marcelo Mollinari, \email{mmollin@ncsu.edu}
#' @export
#'
states_to_visit_mp1 <- function(dat, ph, is.log = TRUE){
  ngam <- choose(dat$ploidy, dat$ploidy/2)
  S <- data.frame(expand.grid(0:(ngam-1),
                              0:(ngam-1))[,2:1],
                  as.numeric(1/(ngam * ngam)))
  Y <- vector("list", length(dat$geno.dose))
  cte<-1
  for(i in 1:nrow(dat$geno.dose))
  {
    a1<-apply(combn(ph[[1]][i,], ncol(ph[[1]])/2), 2, sum)
    a2<-apply(combn(ph[[2]][i,], ncol(ph[[2]])/2), 2, sum)
    st<-apply(expand.grid(a2, a1)[2:1], 1, sum)
    for(j in 1:ncol(dat$geno.dose))
    {
      Y[[cte]] <- data.frame(ind = j,
                             mrk = i,
                             pl.1 = dat$ploidy,
                             pl.2 = dat$ploidy,
                             pl.id = 1,
                             st.p1 = S[st == dat$geno.dose[i,j],1],
                             st.p2 = S[st == dat$geno.dose[i,j],2],
                             st.all = which(st == dat$geno.dose[i,j])-1,
                             emit = S[st == dat$geno.dose[i,j],3],
                             row.names = NULL)
      if(is.log) Y[[cte]]$emit <- log(Y[[cte]]$emit)
      cte <- cte + 1
    }
  }
  Y <- Y %>% bind_rows %>% arrange(ind)
  Y$ind <- as_factor(Y$ind)
  Y$mrk <- as_factor(Y$mrk)
  pls <- matrix(c(dat$ploidy, dat$ploidy), nrow = 1, dimnames = list(paste(dat$ploidy, dat$ploidy, sep = "x"), NULL))
  list(hmm.info = Y, err = 0.0, is.log = TRUE, ploidy.cross.id = pls)
}






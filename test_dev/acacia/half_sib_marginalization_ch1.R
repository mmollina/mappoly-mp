require(mappoly)
require(tidyverse)
setwd("~/repos/current_work/mappoly2/test_dev/acacia/")
dat <- readRDS("data_acacia.rda")
pedigree <- readRDS("pedigree_acacia.rda")
dat.fs.ch1 <- readRDS(file="multiple_fulsib_data_acacia_ch1.rds")
ploidy <- 2


#### Mapping ####
mrks <- lapply(dat.fs.ch1, function(x) x$mrk.names)
a <- Reduce(intersect, mrks)
sm <- a[sort(sample(1:10000, size = 100))]
x <- vector("list", length(dat.fs.ch1))
for(i in 1:length(dat.fs.ch1)){
  cat("full-sub ", names(dat.fs.ch1)[i], "--->", round(i*100/length(dat.fs.ch1),1) ,"% \n")
  d1 <- sample_data(dat.fs.ch1[[i]], type = "marker", selected.mrk = sm)
  s <- make_seq_mappoly(d1, "all")
  tpt <- est_pairwise_rf(s, ll = T, verbose = FALSE)
  v <- tpt$pairwise[names(unlist(lapply(tpt$pairwise, nrow)))]
  w <- lapply(v, function(x) {x[2,,drop = FALSE]})
  x[[i]] <- data.frame(reshape2::melt(w, value.name = "l"), pop = names(dat.fs.ch1)[i])
}
z <- reduce(x, full_join)
z <- cbind(z, str_split_fixed(z$Var2, "-", 2), str_split_fixed(z$pop, "-", 2))[,c(8,9,4,6,7,3)]
colnames(z) <- c("P1", "P2", "mrk_pair", "phaseP1", "phaseP2", "loglike")
w <- reshape2::melt(z, id = c("mrk_pair", "phaseP1", "phaseP2", "loglike"), value.name = "Parent")
w <- w %>%
  mutate(phase = case_when(variable == 'P1' ~ as.numeric(phaseP1),
                           variable == 'P2' ~ as.numeric(phaseP2))) %>%
  select(-one_of(c("phaseP1", "phaseP2", "variable"))) %>%
  group_by(mrk_pair, Parent, phase) %>%
  summarise(sum_loglike = sum(loglike)) %>%
  arrange(Parent) %>%
  arrange(sum_loglike, .by_group = TRUE)%>%
  #filter(length(mrk_pair) > 1)  %>%
  group_by(Parent, mrk_pair) %>%
  mutate(LOD = sum_loglike - lag(sum_loglike, default = sum_loglike[1]))
w

#### Two-points ####
my_func <- function(x, parent, ploidy){
  mrk <- names(x)[2:3]
  y1 <- c(rep(1, x[2]), rep(0, (ploidy-x[2])))
  y2 <- c(rep(1, x[3]), rep(0, (ploidy-x[3])))
  if(length(unique(y1)) == 1 | length(unique(y2)) == 1){
    y <- as.matrix(rbind(y1, y2))
    rownames(y) <- mrk
    colnames(y) <- paste(parent, 1:ncol(y), sep = "_")
    return(y)
  } else {
    v <- perm_tot(y2)
    u <- apply(v, 1, function(x) sum(x + y1 == 2))
    id<-match(unique(u), u)
    y <- as.matrix(rbind(y1, v[id,][u[id]==as.numeric(x[1]),]))
    rownames(y) <- mrk
    colnames(y) <- paste(parent, 1:ncol(y), sep = "_")
    return(y)
  }
}

l <- seq(1,100,length.out = 10)
M <- L <- matrix(NA, 10, 10)
for(i1 in 1:9){
  for(i2 in (i1+1): 10)
  {
   m1 <- l[i1]
   m2 <- l[i2]
   #### Parental phases ####
   y <- w%>% filter(mrk_pair == paste(m1,m2, sep = "-")) %>%
     filter(LOD == 0)
   parent.names <- unique(y$Parent)
   P <- as.data.frame(t(dat[sm[c(m1,m2)], parent.names])) %>%
     rownames_to_column(var = "Parent")
   if(nrow(P) < 2)
     next()
   Z <- left_join(y, P)
   ph <- vector("list", nrow(Z))
   for(i in 1:nrow(Z)){
     ph[[i]] <- my_func(x = Z[i, c(3,6,7)], parent = Z[i, "Parent"], ploidy = 2)
   }
   names(ph) <- parent.names
   ph
   #### Pedigree ####
   all.pop <- apply(expand_grid(parent.names, parent.names), 1, paste0, collapse = "-")
   id.pop <- which(!is.na(match(pedigree$pop, all.pop)))
   pd<-pedigree[id.pop,]
   rownames(pd) <- NULL
   pd <- column_to_rownames(pd, "Ind")
   pd <- data.frame(Par1 = pd[,1], Par2 = pd[,2], pl1 = 2, pl2 = 2, row.names = rownames(pd))
   head(pd)
   pd
   #### Offspring ####
   my_func2 <- function(x, ploidy){
     y <- matrix(0,length(x),ploidy)
     rownames(y) <- names(x)
     for(i in 1:length(x))
       if(x[i]!=0) y[i,1:x[i]] <- 1
     return(y)
   }
   O <- as.data.frame(t(dat[sm[c(m1,m2)], rownames(pd)]))
   of <- apply(O, 1, my_func2, ploidy = 2, simplify = FALSE)
   dat.cur <-  structure(list(offspring = of,
                              phases = ph,
                              pedigree = pd,
                              map = NA,
                              joint.info  = NA),
                         class = "mappoly2.data")
   states <- states_to_visit(dat.cur)
   #save.image("test.Rdata")
   restemp <- hmm_map_reconstruction(ploidy1 = states$ploidy$pl1,
                                     ploidy2 = states$ploidy$pl2,
                                     n.mrk = states$n.mrk,
                                     n.ind = states$n.ind,
                                     haplo = states$states,
                                     emit = states$emit,
                                     rf_vec = rep(0.01, states$n.mrk-1),
                                     verbose = FALSE,
                                     use_H0 = FALSE,
                                     tol = 1e-3)
   restemp.h0 <- hmm_map_reconstruction(ploidy1 = states$ploidy$pl1,
                                        ploidy2 = states$ploidy$pl2,
                                        n.mrk = states$n.mrk,
                                        n.ind = states$n.ind,
                                        haplo = states$states,
                                        emit = states$emit,
                                        rf_vec = rep(0.5, states$n.mrk-1),
                                        verbose = FALSE,
                                        use_H0 = TRUE,
                                        tol = 1e-3)
   #### Results ####
   M[i1, i2] <- M[i2, i1] <- restemp[[2]]
   L[i1, i2] <- L[i2, i1] <- restemp[[1]] - restemp.h0[[1]]
  }
}
heatmap(L, Rowv = NA, Colv = NA)
heatmap(1/log(M), Rowv = NA, Colv = NA)








d1 <- sample_data(dat.fs.ch1[[68]], type = "marker", selected.mrk = sm)
s <- make_seq_mappoly(d1, "all")
tpt <- est_pairwise_rf(s, verbose = FALSE)
m <- rf_list_to_matrix(tpt)
plot(m)
k1 <- 1:5
k2 <- 20:30
round(m$rec.mat[k1, k2], 3)
tpt$pairwise$`1-24`

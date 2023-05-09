require(mappoly)
require(tidyverse)
setwd("~/repos/current_work/mappolymp/test_dev/acacia/")
dat <- readRDS("data_acacia.rda")
pedigree <- readRDS("pedigree_acacia.rda")
dat.fs <- readRDS(file="multiple_fulsib_data_acacia.rda")
ploidy <- 2
#### Selecting half-sib families ####
head(pedigree)
parent <- "AC332"
## selecting individuals with crosses with the same specified parent (half-sib)
ped.hs <- pedigree %>% filter(P1 == parent | P2 == parent) %>% arrange(pop)
n.hs <- ped.hs %>% count(pop)
## selecting half-sib data removing redundant data
dat.hs <- unique(dat[,ped.hs$Ind, drop = FALSE]) ## 73.3k x 74
selcol <- c("snp_name", "ch", "pos", "ref", "alt",  unique(c(ped.hs$P1, ped.hs$P2)))
selrow <- rownames(dat.hs)
dat.hs <- cbind(dat[selrow, selcol], dat.hs)
dat.hs[1:10, 1:20]

## selecting a number of markers (for speed purposes)
#n <- sort(sample(ncol(dat.hs), 500))
n <- 1:1000
image(t(dat.hs[n,ped.hs$Ind]), axes = FALSE)
d <- 1/(nrow(ped.hs)-1)
s <- seq(-d/2, 1+d/2, d)
abline(v = s[cumsum(c(1, n.hs$n))], col = "blue")
#### MAPpoly objects for selected half sib ####
hs <- split(ped.hs, f = ped.hs$pop)
## Removing pops with unknown parents
hs <- hs[str_detect(names(hs), "NA", TRUE)]
## Assembling MAPpoly datasets
dat.hs.mp <- vector("list", length(hs))
npop.init <- length(dat.hs.mp)
names(dat.hs.mp) <- names(hs)
for(i in 1:npop.init){
  cat("pop", i, "---> ")
  ## Selection a full-sib
  dt <- dat.hs[,hs[[i]][,1], drop = FALSE]
  cat(nrow(dt), "markers /", ncol(dt), "individuals")
  dt <- tryCatch(data.frame(snp_name = rownames(dt),
                            P1 = dat.hs[rownames(dt), hs[[i]][1,2], drop = FALSE],
                            P2 = dat.hs[rownames(dt), hs[[i]][1,3], drop = FALSE],
                            ch = dat.hs[rownames(dt), "ch", drop = FALSE],
                            pos = dat.hs[rownames(dt), "pos", drop = FALSE],
                            dt), error = function(e) NA)
  if(all(is.na(dt)))
  {
    cat("\n")
    next()
  }
  dat.hs.mp[[i]] <- tryCatch(table_to_mappoly(dat = dt,
                                              ploidy,
                                              verbose = FALSE,
                                              elim.redundant = FALSE),
                             error = function(e) NA)
  if(!all(is.na(dat.hs[[i]]))){
    dat.hs.mp[[i]]$seq.ref <- dat.hs[rownames(dt), "ref"]
    dat.hs.mp[[i]]$seq.alt <- dat.hs[rownames(dt), "alt"]
  }
  cat("\n")
}
dat.hs.mp <- dat.hs.mp[sapply(dat.hs.mp, function(x) !all(is.na(x)))]
## Proportion of remaining crosses
round(100*length(dat.hs.mp)/npop.init, 1)
saveRDS(dat.hs.mp, file=paste0("half_sibs_", parent,"_data_acacia.rda"))
dat.hs.mp <- readRDS("~/repos/current_work/mappolymp/test_dev/acacia/half_sibs_AC332_data_acacia.rds")





dat.fs
mrks <- lapply(dat.fs, function(x) x$mrk.names)
#### Mapping ####
mrks <- lapply(dat.hs.mp, function(x) x$mrk.names)
a <- Reduce(intersect, mrks)
x <- vector("list", length(dat.hs.mp))
for(i in 1:length(dat.hs.mp)){
  d1 <- sample_data(dat.hs.mp[[i]], type = "marker", selected.mrk = a[1:50])
  p <-match(parent, str_split_fixed(names(dat.hs.mp)[i], "-", 2))
  s <- make_seq_mappoly(d1, "all")
  tpt <- est_pairwise_rf(s, ll = T, verbose = FALSE)
  w <- lapply(tpt$pairwise, function(x) x[2,,drop = FALSE])
  x[[i]] <- data.frame(reshape2::melt(w, value.name = "l"), pop = names(dat.hs.mp)[i])
}
z <- reduce(x, full_join)
z <- cbind(z, str_split_fixed(z$Var2, "-", 2), str_split_fixed(z$pop, "-", 2))[,c(8,9,4,6,7,3)]
colnames(z) <- c("P1", "P2", "mrk_pair", "phaseP1", "phaseP2", "loglike")
head(z)
y <- z %>%
  group_by(mrk_pair, P1, phaseP1) %>%
  summarise(sum.l = sum(loglike))
y


print(y, n = 106)


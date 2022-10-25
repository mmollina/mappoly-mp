# ## Saving Pedigree
# load("~/repos/acacia/Pedigree_G.reduced.RData")
# pedigree <- Pedigree_G.reduced
# colnames(pedigree) <- c("Ind", "P1", "P2")
# pedigree$Ind <- paste0("Ind_", pedigree$Ind)
#
# saveRDS(pedigree, file = "test_dev/pedigree_acacia.rda")
# 
# ## Saving Data
# load("~/repos/acacia/MAP.non.redundant.Rdata")
# load("~/repos/acacia/dataHD.non.redundant.Rdata")
# info <- map.non.redundant[,c(1,2,3,4,5)]
# rownames(info) <- info[,1]
# dat <- t(dataHD.non.redundant)
# dat <- unique(dat)
# dat <- data.frame(snp_name = rownames(dat),
#                   ch = info[rownames(dat),2],
#                   pos = info[rownames(dat),3],
#                   ref = info[rownames(dat),4],
#                   alt = info[rownames(dat),5],
#                   dat)
# colnames(dat) <- str_replace_all(string = colnames(dat), pattern = "X", replacement = "Ind_")
# saveRDS(dat, file="test_dev/data_acacia.rda")
require(mappoly)
require(stringr)
require(ggplot2)
require(dplyr)
dat <- readRDS("test_dev/data_acacia.rda")
pedigree <- readRDS("test_dev/pedigree_acacia.rda")

#### Converting to MAPpoly ####
ploidy <-2
## labeling populations
pedigree$pop <- apply(pedigree[, 2:3], 1, paste0, collapse="-")
## Split into full-sib
fs <- split(pedigree, f = pedigree$pop)
## Removing pops with unknown parents
fs <- fs[str_detect(names(fs), "NA", TRUE)]
## Assembling MAPpoly datasets
dat.fs <- vector("list", length(fs))
npop.init <- length(dat.fs)
names(dat.fs) <- names(fs)
for(i in 1:npop.init){
  cat("pop", i, "---> ")
  ## Selection a full-sib
  dt <- dat[,fs[[i]][,1], drop = FALSE]
  dt <- unique(dt)
  cat(nrow(dt), "markers /", ncol(dt), "individuals")
  dt <- tryCatch(data.frame(snp_name = rownames(dt),
                            P1 = dat[rownames(dt), fs[[i]][1,2], drop = FALSE],
                            P2 = dat[rownames(dt), fs[[i]][1,3], drop = FALSE],
                            ch = dat[rownames(dt), "ch", drop = FALSE],
                            pos = dat[rownames(dt), "pos", drop = FALSE],
                            dt), error = function(e) NA)
  if(all(is.na(dt)))
  {
    cat("\n")
    next()    
  }
  dat.fs[[i]] <- tryCatch(table_to_mappoly(dt, ploidy, verbose = FALSE), error = function(e) NA)
  if(!all(is.na(dat.fs[[i]]))){
    dat.fs[[i]]$seq.ref <- dat[rownames(dt), "ref"]
    dat.fs[[i]]$seq.alt <- dat[rownames(dt), "alt"]
  }
  cat("\n")
}
dat.fs <- dat.fs[sapply(dat.fs, function(x) !all(is.na(x)))]
## Proportion of remaining crosses 
round(100*length(dat.fs)/npop.init, 1)

#### Plot ####
rcross <- str_split_fixed(names(dat.fs), "-", 2)
rcross <- expand.grid(unique(rcross[,1]), unique(rcross[,2]))
dimnames(rcross) <- list(apply(rcross, 1, paste0, collapse = "-"), c("P1", "P2"))
rcross$n.ind <- NA
rcross$n.mrk <- NA
n.ind <- sapply(dat.fs, function(x) x$n.ind)
rcross[names(n.ind),"n.ind"] <- n.ind
n.mrk <- sapply(dat.fs, function(x) x$n.mrk)
rcross[names(n.mrk),"n.mrk"] <- n.mrk
P1 <- rcross %>% group_by(P1) %>% summarise(tot.p1 = sum(n.ind, na.rm = TRUE))
P2 <- rcross %>% group_by(P2) %>% summarise(tot.p2 = sum(n.ind, na.rm = TRUE))
head(rcross)
rcross2 <- rbind(data.frame(P1 = P1$P1,
                           P2 = "TotP2",
                           n.ind = P1$tot.p1,
                           n.mrk = NA),
                data.frame(P1 = "TotP1",
                           P2 = P2$P2,
                           n.ind = P2$tot.p2,
                           n.mrk = NA), rcross)

ggplot(rcross2,
       aes(x = str_to_title(P1), 
           y = str_to_title(P2),
           color = n.mrk,
           size = n.ind)) +
  geom_point() +
  geom_text(aes(label = n.ind), 
            colour = "white", 
            size = 3) +
  scale_x_discrete(position = "top") +
  scale_size_continuous(range = c(5, 15)) + # Adjust as required.
  scale_color_gradientn(colors = c("red", "darkblue"), na.value = "darkgray") +
  labs(x = NULL, y = NULL)

#### Selecting half-sib families ####
head(pedigree)
parent <- "AC332"
id <- which(apply(pedigree, 1, function(x) any(x[2:3] == parent)))
a <- pedigree[id,]
dt <- dat[,a[,1], drop = FALSE]
dt <- t(dt)
M<-cor(dt)
dt <- unique(dt)

cat(nrow(dt), "markers /", ncol(dt), "individuals")
dt <- tryCatch(data.frame(snp_name = rownames(dt),
                          P1 = dat[rownames(dt), a[1,2], drop = FALSE],
                          P2 = dat[rownames(dt), a[1,3], drop = FALSE],
                          ch = dat[rownames(dt), "ch", drop = FALSE],
                          pos = dat[rownames(dt), "pos", drop = FALSE],
                          dt), error = function(e) NA)
if(all(is.na(dt)))
{
  cat("\n")
  next()    
}
dat.fs[[i]] <- tryCatch(table_to_mappoly(dt, ploidy, verbose = FALSE), error = function(e) NA)
if(!all(is.na(dat.fs[[i]]))){
  dat.fs[[i]]$seq.ref <- dat[rownames(dt), "ref"]
  dat.fs[[i]]$seq.alt <- dat[rownames(dt), "alt"]
}
cat("\n")


id <- which(str_detect(names(dat.fs), parent))
d1 <- dat.fs[[id[1]]]
s <- make_seq_mappoly(d1, "all")
tpt <- est_pairwise_rf(s, ncpus = 32)
tpt$pairwise$`1-2`


m <- rf_list_to_matrix(tpt, 2,2,.1)
plot(m)
g <- group_mappoly(m, expected.groups = 13, comp.mat = T)
g
s1 <- make_seq_mappoly(g, 1)
tpt1 <- make_pairs_mappoly(tpt, s1)
map <- est_rf_hmm_sequential(input.seq = s1,
                             twopt = tpt1,
                             extend.tail = 20,
                             sub.map.size.diff.limit = 8)
map.up <- est_full_hmm_with_global_error(input.map = map,
                                         error = 0.1,
                                         verbose = TRUE)
print(map.up, detailed = T)
plot(map.up)




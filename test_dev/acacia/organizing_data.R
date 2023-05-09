# ## Saving Pedigree
setwd("~/repos/current_work/mappolymp/test_dev/acacia/")
# load("Pedigree_G.reduced.RData")
# pedigree <- Pedigree_G.reduced
# colnames(pedigree) <- c("Ind", "P1", "P2")
# pedigree$Ind <- paste0("Ind_", pedigree$Ind)
#
# saveRDS(pedigree, file = "test_dev/pedigree_acacia.rda")
#
# ## Saving Data
# load("MAP.non.redundant.Rdata")
# load("dataHD.non.redundant.Rdata")
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
# saveRDS(dat, file="data_acacia.rda")
require(mappoly)
require(stringr)
require(ggplot2)
require(dplyr)
dat <- readRDS("data_acacia.rda")
pedigree <- readRDS("pedigree_acacia.rda")
dat[1:10, 50:60]
pedigree[25:35,]
#### Converting to MAPpoly ####
ploidy <- 2
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
saveRDS(dat.fs, file="multiple_fulsib_data_acacia.rda")

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

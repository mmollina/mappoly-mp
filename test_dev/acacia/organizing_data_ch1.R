# ## Saving Pedigree
setwd("~/repos/current_work/mappoly2/test_dev/acacia/")
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
require(tidyverse)
dat <- readRDS("data_acacia.rda")
pedigree <- readRDS("pedigree_acacia.rda")
#### Converting to MAPpoly ####
ploidy <- 2
## labeling populations
pedigree$pop <- apply(pedigree[, 2:3], 1, paste0, collapse="-")
## Split into full-sib
fs <- split(pedigree, f = pedigree$pop)
## Removing pops with unknown parents
fs <- fs[str_detect(names(fs), "NA", TRUE)]
## Assembling MAPpoly datasets
dat.fs.ch1 <- vector("list", length(fs))
npop.init <- length(dat.fs.ch1)
names(dat.fs.ch1) <- names(fs)
for(i in 1:npop.init){
  cat("pop", i, "---> ")
  ## Selection a full-sib
  dt <- dat[dat$ch == 1,fs[[i]][,1], drop = FALSE]
  cat(nrow(dt), "markers /", ncol(dt), "individuals ---> ")
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
  dat.fs.ch1[[i]] <- tryCatch(mp2_csv_to_mappoly(dt),
                              error = function(e) NA)
  #dat.fs.ch1[[i]] <- tryCatch(mappoly::table_to_mappoly(dt, ploidy = 2,
  #                                                      elim.redundant = FALSE,
  #                                                      verbose = FALSE),
  #                            error = function(e) NA)
  dat.fs.ch1[[i]]$seq.ref <- dat[rownames(dt), "ref"]
  dat.fs.ch1[[i]]$seq.alt <- dat[rownames(dt), "alt"]
  names(dat.fs.ch1[[i]]$seq.alt) <- names(dat.fs.ch1[[i]]$seq.ref) <- dat.fs.ch1[[1]]$mrk.names
  cat(dat.fs.ch1[[i]]$n.mrk)
  cat("\n")
}
dat.fs.ch1 <- dat.fs.ch1[sapply(dat.fs.ch1, function(x) !all(is.na(x)))]
## Proportion of remaining crosses
round(100*length(dat.fs.ch1)/npop.init, 1)
saveRDS(dat.fs.ch1, file="multiple_fulsib_data_acacia_ch1.rds")

require(mappoly2)
require(mappoly)
require(tidyverse)
ploidy <- ploidy1 <- 4
ploidy2 <- 4
n.mrk <- 10
n.ind <- 100
cm.map <- seq(0, 10, length.out = n.mrk)

#### Simulating ####
h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
                                           n.mrk = n.mrk,
                                           alleles = 0:1,
                                           lambda = 1,
                                           shuffle.homolog = TRUE)
op <- par(mfrow = c(1,2))
image(h1, axes = FALSE, xlab = "Markers", ylab = "homologs")
abline(h = seq(-1/(2*ploidy1 - 2), 1 + 1/(2*ploidy1- 2),
               length.out = ploidy1+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h1))
axis(2, at = seq(0,1, length.out = ploidy1), labels = colnames(h1))

h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
                                           n.mrk = n.mrk,
                                           alleles = 0:1,
                                           lambda = 1,
                                           shuffle.homolog = TRUE)
image(h2, axes = FALSE, xlab = "Markers", ylab = "homologs")
abline(h = seq(-1/(2*ploidy2 - 2), 1 + 1/(2*ploidy2 - 2),
               length.out = ploidy2+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h2))
axis(2, at = seq(0,1, length.out = ploidy2), labels = colnames(h2))
par(op)
f1 <- simulate_cross(h1, h2, cm.map)

#### MAPPoly ####
dose.geno <- apply(f1$offspring, c(1,3), sum)
colnames(dose.geno) <- paste0("Ind", 1:ncol(dose.geno))
dat <- data.frame(snp.name = paste0("M", 1:n.mrk),
                  P1 = apply(f1$ph1 , 1, sum),
                  P2 = apply(f1$ph2 , 1, sum),
                  sequence = 1,
                  sequence.pos = 1:n.mrk,
                  dose.geno)

dat.mp <- table_to_mappoly(dat, ploidy)
s <- make_seq_mappoly(dat.mp, "all")
tpt <- est_pairwise_rf(s)
map <- est_rf_hmm_sequential(s,tpt, tol.final = 10e-5)
print(map, detailed = TRUE)
plot(map)
plot_compare_haplotypes(ploidy,
                        ph_matrix_to_list(f1$ph1[map$info$mrk.names,]),
                        ph_matrix_to_list(f1$ph2[map$info$mrk.names,]),
                        map$maps[[1]]$seq.ph$P,
                        map$maps[[1]]$seq.ph$Q)
mp <- round(cumsum(mappoly::imf_h(c(0, map$maps[[1]]$seq.rf))),2)
#### States to visit ####
ngam1 <- choose(ploidy1, ploidy1/2)
ngam2 <- choose(ploidy2, ploidy2/2)
S <- cbind(ploidy1, ploidy2, 1, expand.grid(0:(ngam1-1), 0:(ngam2-1))[,2:1],
           0:((ngam1*ngam2)-1), 0)
colnames(S) <- c("pl1", "pl2", "pl.id", "st.p1", "st.p2", "st.all", "emit")
Y <- NULL
info.mrk <- as.numeric(str_remove_all(map$info$mrk.names, "M"))
for(i in 1:dim(f1$offspring)[3]){
  for(j in info.mrk){
    A1 <- combn(as.numeric(f1$ph1[j,]), ploidy1/2)
    A2 <- combn(as.numeric(f1$ph2[j,]), ploidy2/2)
    w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                   apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
    a1 <- apply(str_split_fixed(w,"_", (ploidy1 + ploidy2)/2),1,sort)
    id <- which(apply(a1, 2, function(x) all(x == sort(f1$offspring[j,,i]))))
    tmp <- S[id,]
    tmp$emit[] <- log(1/length(tmp$emit))
    Y<-rbind(Y, cbind(ind = paste0("Ind_", i), mrk = paste0("M_", j), tmp))
  }
}
Y$ind <- as.character(Y$ind)
Y$mrk <- as.character(Y$mrk)
sts <- list(hmm.info = Y,
            err = 0,
            is.log = TRUE,
            ploidy.cross.id = matrix(c(ploidy1,ploidy2),1,2, dimnames = list(paste(ploidy1, ploidy2, sep = "x"), NULL)))

#### MAPpoly haplotype ####
mrks <- as.character(unique(sts$hmm.info$mrk))
emit <- h <- vector("list", length(mrks))
names(emit)  <- names(h) <- mrks
for(i in mrks)
{
  st.temp <- sts$hmm.info %>%
    filter(mrk == i)
  ind <- as.character(unique(st.temp$ind))
  emit.temp <- htemp <- vector("list", length(ind))
  names(emit.temp) <- names(htemp) <- ind
  for(j in ind){
    htemp[[j]] <- as.matrix(st.temp %>%
                              filter(ind == j) %>%
                              select("st.p1", "st.p2"))
    emit.temp[[j]] <- exp(as.matrix(st.temp %>%
                              filter(ind == j) %>%
                              select("emit")))
  }
  h[[i]] <- htemp
  emit[[i]] <- emit.temp
}
restemp <- mappoly:::est_haplo_hmm(ploidy = ploidy,
                                   n.mrk = length(h),
                                   n.ind = length(h[[1]]),
                                   haplo = h,
                                   emit = emit,
                                   rf_vec = rep(0.01, length(h)-1),
                                   verbose = FALSE,
                                   use_H0 = FALSE,
                                   tol = 10e-5)
mp.hap <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
#### Comparing 0####
Y0 <- rbind(cm.map[info.mrk],mp,mp.hap)
rownames(Y0) <- c("simulation", "mappoly", "mappoly.haplotype")
Y0
#### Checking the effect of changing emission function ####
emit2 <- emit
for(i in mrks)
{
  for(j in ind){
    emit2[[i]][[j]] <- emit[[i]][[j]] + runif(length(emit[[i]][[j]]), -0.01, 0.01)
  }
}
restemp <- mappoly:::est_haplo_hmm(ploidy = ploidy,
                                   n.mrk = length(h),
                                   n.ind = length(h[[1]]),
                                   haplo = h,
                                   emit = emit2,
                                   rf_vec = rep(0.01, length(h)-1),
                                   verbose = FALSE,
                                   use_H0 = FALSE,
                                   tol = 10e-5)
mp.hap.e <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
#### Comparing 1####
Y1 <- rbind(cm.map[info.mrk],mp,mp.hap,mp.hap.e)
rownames(Y1) <- c("simulation", "mappoly", "mappoly.haplotype", "mappoly.perturbed.emssion")
Y1

#### MAPpoly2 - same C++ code as mappoly####
restemp <- hmm_map_reconstruction(ploidy = ploidy,
                                   n.mrk = length(h),
                                   n.ind = length(h[[1]]),
                                   haplo = h,
                                   emit = emit,
                                   rf_vec = rep(0.01, length(h)-1),
                                   verbose = FALSE,
                                   use_H0 = FALSE,
                                   tol = 10e-5)
mp2 <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
#### Comparing 2####
Y2 <- rbind(cm.map[info.mrk],mp,mp.hap,mp.hap.e, mp2)
rownames(Y2) <- c("simulation", "mappoly", "mappoly.haplotype", "mappoly.perturbed.emssion", "mappoly2")
Y2
##############Multiallelism################
ploidy1 <- 4
ploidy2 <- 4
n.mrk <- 10
n.ind <- 200
al1 <- 1:4
al2 <- 5:8
cm.map <- seq(0, 5, length.out = n.mrk)
####Simulation ####
h1 <- simulate_multiallelic_homology_group(ploidy = ploidy1,
                                           n.mrk = n.mrk,
                                           alleles = al1,
                                           shuffle.homolog = TRUE)
h2 <- simulate_multiallelic_homology_group(ploidy = ploidy2,
                                           n.mrk = n.mrk,
                                           alleles = al2,
                                           shuffle.homolog = TRUE)
op <- par(mfrow = c(1,2))
col1 <- sort(unique(as.vector(h1)))
col2 <- sort(unique(as.vector(h2)))
pal <- RColorBrewer::brewer.pal(length(unique(c(col1, col2))), "Set1")
col1 <- pal[col1]
col2 <- pal[col2]
fields::image.plot(h1, axes = FALSE, xlab = "Markers", ylab = "homologs", col = col1)
abline(h = seq(-1/(2*ploidy1 - 2), 1 + 1/(2*ploidy1- 2),
               length.out = ploidy1+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h1))
axis(2, at = seq(0,1, length.out = ploidy1), labels = colnames(h1))
fields::image.plot(h2, axes = FALSE, xlab = "Markers", ylab = "homologs", col = col2)
abline(h = seq(-1/(2*ploidy2 - 2), 1 + 1/(2*ploidy2 - 2),
               length.out = ploidy2+1),
       lwd = 2)
axis(1, at = seq(0,1, length.out = n.mrk), labels = rownames(h2))
axis(2, at = seq(0,1, length.out = ploidy2), labels = colnames(h2))
par(op)
f1 <- simulate_cross(h1, h2, cm.map)

####States to visit####
ngam1 <- choose(ploidy1, ploidy1/2)
ngam2 <- choose(ploidy2, ploidy2/2)
S <- cbind(ploidy1, ploidy2, 1,
           expand.grid(0:(ngam1-1), 0:(ngam2-1))[,2:1],
           0:((ngam1*ngam2)-1), 0)
colnames(S) <- c("pl1", "pl2", "pl.id", "st.p1", "st.p2", "st.all", "emit")
Y <- NULL
b1 <- apply(f1$ph1, 1, function(x) length(unique(x)))
b2 <- apply(f1$ph2, 1, function(x) length(unique(x)))
info.mrk <- as.numeric(str_remove_all(names(which(!(b1==1 & b2==1))), "M"))
for(i in 1:dim(f1$offspring)[3]){
  for(j in info.mrk){
    A1 <- combn(as.numeric(f1$ph1[j,]), ploidy1/2) ## Ordered vector ---> phased
    A2 <- combn(as.numeric(f1$ph2[j,]), ploidy2/2) ## Ordered vector ---> phased
    w <- kronecker(apply(A1, 2, paste0, collapse = "_"),
                   apply(A2, 2, paste0, collapse = "_"), paste, sep = "_")
    a1 <- apply(str_split_fixed(w,"_", (ploidy1 + ploidy2)/2),1,sort)
    id <- which(apply(a1, 2, function(x) all(x == sort(f1$offspring[j,,i]))))
    tmp <- S[id,]
    tmp$emit[] <- log(1/length(tmp$emit))
    Y<-rbind(Y, cbind(ind = paste0("Ind_", i), mrk = paste0("M_", j), tmp))
  }
}
sts <- list(hmm.info = Y,
            err = 0,
            is.log = TRUE,
            ploidy.cross.id = matrix(c(ploidy1,ploidy2),1,2, dimnames = list(paste(ploidy1, ploidy2, sep = "x"), NULL)))
mrks <- as.character(unique(sts$hmm.info$mrk))
emit <- h <- vector("list", length(mrks))
names(emit)  <- names(h) <- mrks
for(i in mrks)
{
  st.temp <- sts$hmm.info %>%
    filter(mrk == i)
  ind <- as.character(unique(st.temp$ind))
  emit.temp <- htemp <- vector("list", length(ind))
  names(emit.temp) <- names(htemp) <- ind
  for(j in ind){
    htemp[[j]] <- as.matrix(st.temp %>%
                              filter(ind == j) %>%
                              select("st.p1", "st.p2"))
    emit.temp[[j]] <- exp(as.matrix(st.temp %>%
                                      filter(ind == j) %>%
                                      select("emit")))
  }
  h[[i]] <- htemp
  emit[[i]] <- emit.temp
}

#### MAPpoly2 - same C++ code as mappoly####
restemp <- hmm_map_reconstruction(ploidy = ploidy,
                                  n.mrk = length(h),
                                  n.ind = length(h[[1]]),
                                  haplo = h,
                                  emit = emit,
                                  rf_vec = rep(0.01, length(h)-1),
                                  verbose = FALSE,
                                  use_H0 = FALSE,
                                  tol = 10e-5)
mp.multiallelic <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)





#### Firts bias assessment ####
require(mappoly2)
require(mappoly)
require(tidyverse)
source("sand/myfunc.R")
ploidy <- ploidy1 <- 4
ploidy2 <- 4
n.mrk <- 10
n.ind <- 200
al1 <- 1:4
al2 <- 5:8
cm.map <- seq(0, 5, length.out = n.mrk)

dev.off()
plot(0, xlim = c(0,5), ylim = c(0,5), type = "n")
abline(a = 0, 1, col = 2, lwd = 3)
for(ble in 1:4){
  cat("simulation: ", ble, "\n")
  z<-myfunc()
  a <- lm(z[[2]]~z[[1]])
  points(z[[2]]~z[[1]], cex = .5, col = "lightgray")
  abline(a, lwd = .5, col = 4)
}










require(mappoly2)
require(mappoly)
require(tidyverse)
p <- c(4,4)
n.ind <- 200
n.mrk <- 3
map.length <- 5
#### Simulating data using mappoly2 ####
names(p) <- c("P1", "P2")
cm <- matrix(c("P1","P2"),
             ncol = 2, byrow = T)
flag <- TRUE
while(flag){
  D <- simulate_multiple_crosses(ploidy = p, #four parents
                                 cross.mat = cm,
                                 n.ind = rep(n.ind, nrow(cm)), # per cross
                                 n.mrk= rep(n.mrk, length(p)), # per parent
                                 alleles = list(c(1:0), c(1:0)),
                                 map.length = map.length, lambda = c(0.5, 0.5))

  flag<- !(all(apply(D$phases$P1 , 1, sum) > 0) &
           all(apply(D$phases$P2 , 1, sum) == 0))
  #flag <- !(all(D$phases$P1 == matrix(c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), 4, 4)) &
  #          all(D$phases$P2 == matrix(c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), 4, 4)))
  flag <- FALSE
}
dev.off()
par(mfrow = c(2,2))
image(t(D$phases$P1), axes = F, xlab = "homologs", ylab = "markers", main = "P1")
image(t(D$phases$P2), axes = F, xlab = "homologs", ylab = "markers", main = "P2")
plot(x = D$map, y = rep(1, n.mrk),  xlim = c(0, map.length), pch = 15, col = 2)
abline(h = 1)
#### Converting and analyzing data using mappoly ####
bipar.dat <- mappoly2_to_mappoly(D)
#plot(bipar.dat[[1]])
dat <- bipar.dat[[1]]
s <- make_seq_mappoly(dat, "all")
tpt <- est_pairwise_rf(s)
map <- est_rf_hmm_sequential(s,tpt, tol.final = 10e-5)
print(map, detailed = TRUE)
mp <- round(cumsum(mappoly::imf_h(c(0, map$maps[[1]]$seq.rf))),2)

#### Using mappoly2 ####
states.hmm <- states_to_visit(input.data = D, err = 0.0, is.log = TRUE)
hist(states.hmm$hmm.info$emit)
x1 <- est_map_R(states.hmm,tol = 10e-5, verbose = F)
mp2 <- round(cumsum(mappoly::imf_h(c(0, x1[[2]]))), 2)
#### Comparing 1 ####
Y1 <- rbind(D$map,mp,mp2)
rownames(Y1) <- c("simulation", "mappoly", "mappoly2")
Y1

#### Using mappoly's 'est_haplo_hmm' ####
hl <- states_to_visit(input.data = D, err = 0.0, is.log = TRUE)
mrks <- as.character(unique(hl$hmm.info$mrk))
h <- vector("list", length(mrks))
names(h) <- mrks
for(i in mrks)
{
  st.temp <- hl$hmm.info %>%
    filter(mrk == i)
  ind <- as.character(unique(st.temp$ind))
  htemp <- vector("list", length(ind))
  names(htemp) <- ind
  for(j in ind)
    htemp[[j]] <- as.matrix(st.temp %>%
                              filter(ind == j) %>%
                              select("st.p1", "st.p2"))
  h[[i]] <- htemp
}
restemp <- mappoly:::est_haplo_hmm(ploidy = p[1],
                                   n.mrk = length(h),
                                   n.ind = length(h[[1]]),
                                   haplo = h,
                                   rf_vec = rep(0.01, length(h)-1),
                                   verbose = FALSE,
                                   use_H0 = FALSE,
                                   tol = 10e-5)
mp.hap <- round(cumsum(mappoly::imf_h(c(0, restemp[[2]]))), 2)
#### Comparing 2####
Y2 <- rbind(D$ma,mp,mp.hap,mp2)
rownames(Y2) <- c("simulation", "mappoly", "mappoly.haplotype", "mappoly2")
Y2





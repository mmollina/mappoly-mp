require(mappolymp)
require(mappoly)
#### Functions ####
build_map <- function(input.data,
                      is.biallelic = TRUE,
                      verbose = TRUE,
                      use_H0 = FALSE,
                      tol = 10e-4){
  pedigree <- input.data$pedigree
  a <- apply(pedigree[,1:2], 1, paste0, collapse = "x")
  id <- 1:length(unique(a))
  names(id) <- unique(a)
  pedigree$id <- id[a]
  PH <- input.data$phases
  for(i in 1:length(PH))
    PH[[i]]<-as.matrix(PH[[i]])
  x <- 1:length(PH)
  names(x) <- names(PH)
  pedigree$Par1 <- x[pedigree$Par1]
  pedigree$Par2 <- x[pedigree$Par2]
  pedigree <- as.matrix(pedigree)
  GENO <- input.data$offspring
  for(i in 1:length(GENO))
    GENO[[i]]<-as.matrix(GENO[[i]])
  if(is.biallelic){
    G<-sapply(GENO, function(x) apply(x, 1, sum))
    w <- mappolymp:::vs_biallelic_Rcpp(PH, G, pedigree)
  } else {
    w <- mappolymp:::vs_multiallelic_Rcpp(PH, GENO, pedigree)
  }
  ploidy1_id <- input.data$pedigree$pl1/2-1
  ploidy2_id <- input.data$pedigree$pl2/2-1
  n.mrk <- nrow(input.data$phases[[1]])
  n.ind <- length(input.data$offspring)
  res.temp <-
    .Call("est_hmm_map",
          ploidy1_id,
          ploidy2_id,
          n.mrk,
          n.ind,
          w$states,
          w$emit,
          rep(0.01, n.mrk-1),
          verbose,
          tol,
          use_H0,
          PACKAGE = "mappolymp")
  structure(list(n.mrk = n.mrk,
                 n.ind = n.ind,
                 rf = res.temp[[2]],
                 loglike = res.temp[[1]],
                 states = w$states,
                 emit = w$emit,
                 #phases = input.data$merged.phases,
                 pedigree = input.data$pedigree), class = "multi.mappoly.map")
}
#### Simulation ####
ploidy.vec <- c(4, 2, 6, 4) #three parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4")
cross.mat <- matrix(c("P1","P2",
                      "P3","P1",
                      "P4","P2",
                      "P3","P4"),
                    ncol = 2, byrow = T)
n.mrk <- c(100,100,100,100) #per parent
map.length <- 50 #in centimorgans
alleles <- list(P1 = c(0:1),
                P2 = c(0:1),
                P3 = c(0:1),
                P4 = c(0:1))
n.ind <- c(100, 100, 100, 100) #per cross
input.data <- simulate_multiple_crosses(ploidy.vec,
                                        cross.mat,
                                        n.ind,
                                        n.mrk,
                                        alleles,
                                        map.length)
input.data
#### old - R based ####
w <- states_to_visit(input.data)
map1 <- hmm_map_reconstruction(w, tol = 10e-4, verbose = TRUE)
plot_map_list_consensus(map.list = NULL, map1)
#### new - C++ based ####
map2 <- build_map(input.data, tol = 10e-4)
plot_map_list_consensus(map.list = NULL, map2)
all(map1$rf - map2$rf < 0.0001)





#### Development section ####
my_Rfunc1 <- function(upd){
  L <- vector("list", nrow(upd))
  for(i in 1:nrow(upd)){
    ngam1 <- choose(upd[i,3], upd[i,3]/2)
    ngam2 <- choose(upd[i,4], upd[i,4]/2)
    S <- expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1]
    L[[i]] <- S
  }
  return(L)
}
my_Rfunc2 <- function(v1 = c(0,0,1,2),
                      v2 = c(0,1)){
  x1 <- combn(v1, length(v1)/2)
  x2 <- combn(v2, length(v2)/2)
  m <- matrix(0, ncol(x1)*ncol(x2), nrow(x1) + nrow(x2))
  for(i in 1:ncol(x1)){
    for(j in 1:ncol(x2)){
      m[(i-1) * ncol(x2) + j,] <- sort(c(x1[,i], x2[,j]))
    }
  }
  return(m)
}



#Convert the following R function into Rcpp:
vs <- function(PH, ## list of parental phases
               GENO, ## list of offspring genotypes
               unique_pop_mat, ## data frame containing the unique populations
               pedigree) ## data frame containing the pedigree information
{
  L <- vector("list", nrow(unique_pop_mat))
  for(i in 1:nrow(unique_pop_mat)){
    ngam1 <- choose(unique_pop_mat[i,3], unique_pop_mat[i,3]/2)
    ngam2 <- choose(unique_pop_mat[i,4], unique_pop_mat[i,4]/2)
    S <- expand.grid(0:(ngam2-1), 0:(ngam1-1))[,2:1]
    L[[i]] <- S
  }
  mrk_names <- rownames(PH[[1]])
  H <- vector("list", length(mrk_names))
  for(i in 1:length(H))
    H[[i]] <- vector("list", nrow(pedigree))

  for(pop_id in 1:nrow(unique_pop_mat)){


    ind_id <- which(pedigree[,"id"] == pop_id)


    for(k in 1:length(mrk_names)){

      v1 <- as.numeric(PH[[unique_pop_mat[pop_id == unique_pop_mat[,"id"],1]]][k,])
      v2 <- as.numeric(PH[[unique_pop_mat[pop_id == unique_pop_mat[,"id"],2]]][k,])


      if(any(is.na(v1)) | any(is.na(v2))){
        for(j in ind_id){
          H[[k]][[j]] <- L[[pop_id]]
        }
      } else{
        x1 <- combn(v1, length(v1)/2)
        x2 <- combn(v2, length(v2)/2)
        x <- matrix(0, ncol(x1)*ncol(x2), nrow(x1) + nrow(x2))
        for(i in 1:ncol(x1)){
          for(j in 1:ncol(x2)){
            x[(i-1) * ncol(x2) + j,] <- sort(c(x1[,i], x2[,j]))
          }
        }
        for(j in ind_id){
          a <- GENO[[j]][k,]
          if(any(is.na(a))){
            H[[k]][[j]] <- L[[pop_id]]
          } else {
            g <- sort(a)
            y <- NULL
            for(i in 1:nrow(x)){
              if(all(x[i,] == g))
                y <- c(y, i)
            }
            H[[k]][[j]] <- L[[pop_id]][y,]
          }
        }
      }
    }
  }
  return(H)
}

require(mappolymp)
system.time(w1 <- vs(PH, GENO, unique_pop_mat, pedigree))
w2 <- states_to_visit(input.data)
system.time(w3 <- mappolymp:::vs_multiallelic_Rcpp(PH, GENO, pedigree))

x<-0
for(i in 1:length(w1))
  for(j in 1:length(w1[[i]]))
    x <- x + sum(w2$states[[i]][[j]] - w3$states[[i]][[j]])
x

### Bi-allelic ###
G<-sapply(GENO, function(x) apply(x, 1, sum))
G[1:5, 1:5]
system.time(w4 <- mappolymp:::vs_biallelic_Rcpp(PH, G, pedigree))
x<-0
for(i in 1:length(w1))
  for(j in 1:length(w1[[i]]))
    x <- x + sum(w2$states[[i]][[j]] - w4$states[[i]][[j]])
x

w4$states[[4]][[1]]
w2$states[[4]][[1]]


w4$emit[[4]][[1]]
w2$emit[[4]][[1]]






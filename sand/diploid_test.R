require(mappoly2)
require(onemap)

my_func<-function(D){
  mrk.ind <- dim(D$dat$P1xP2)[c(1,3)]
  R <- NULL
  for(j in 1:mrk.ind[1]){
    mrk.name <- rownames(D$phases$P1)[j]
    p1 <- unlist(D$phases$P1[j,])
    l1 <- length(unique(p1))
    p2 <- unlist(D$phases$P2[j,])
    l2 <- length(unique(p2))
    if(l1 == 2 & l2 == 2){
      if(all(sort(p1) == sort(p2))){
        type <- "B3.7"
        x<-matrix(c(p1[1], p2[1],
                    p1[1], p2[2],
                    p1[2], p2[1],
                    p1[2], p2[2]), byrow = TRUE, ncol = 2)
        id <- apply(x, 1, paste0, collapse = "")
        y <-c("a", "ab", "ab", "b")
        names(y) <- id
        d <- apply(D$dat$P1xP2[j,,], 2, paste0, collapse = "")
        R <- rbind(R, paste(paste0("*M",mrk.name), type, paste(y[d], collapse = " ")))
      }
      else{
        type <- "A.1"
        x<-matrix(c(p1[1], p2[1],
                    p1[1], p2[2],
                    p1[2], p2[1],
                    p1[2], p2[2]), byrow = TRUE, ncol = 2)
        id <- apply(x, 1, paste0, collapse = "")
        y <-c("ac", "ad", "bc", "bd")
        names(y) <- id
        d <- apply(D$dat$P1xP2[j,,], 2, paste0, collapse = "")
        R <- rbind(R, paste(paste0("*M",mrk.name), type, paste(y[d], collapse = " ")))
      }
    }
    if(l1 == 2 & l2 == 1){
      type <- "D1.10"
      x<-matrix(c(p1[1], p2[1],
                  p1[1], p2[2],
                  p1[2], p2[1],
                  p1[2], p2[2]), byrow = TRUE, ncol = 2)
      id <- apply(x, 1, paste0, collapse = "")
      y <-c("a", "a", "ab", "ab")
      names(y) <- id
      d <- apply(D$dat$P1xP2[j,,], 2, paste0, collapse = "")
      R <- rbind(R, paste(paste0("*M",mrk.name), type, paste(y[d], collapse = " ")))
    }
    if(l1 == 1 & l2 == 2){
      type <- "D2.15"
      x<-matrix(c(p1[1], p2[1],
                  p1[1], p2[2],
                  p1[2], p2[1],
                  p1[2], p2[2]), byrow = TRUE, ncol = 2)
      id <- apply(x, 1, paste0, collapse = "")
      y <-c("a", "ab", "a", "ab")
      names(y) <- id
      d <- apply(D$dat$P1xP2[j,,], 2, paste0, collapse = "")
      R <- rbind(R, paste(paste0("*M",mrk.name), type, paste(y[d], collapse = " ")))
    }
  }
  write("data type outcross", file = "diploid_test.txt")
  write(c(mrk.ind[2], nrow(R),0,0,0), file = "diploid_test.txt", append = TRUE)
  cat(rownames(D$pedigree), "\n", file = "diploid_test.txt", append = TRUE)
  write(R, file = "diploid_test.txt", append = TRUE)
  (dat.onemap<-read_onemap("diploid_test.txt", verbose = FALSE))
  twopt <- rf_2pts(dat.onemap, verbose = FALSE)
  markers <- make_seq(twopt,"all") # correct phases
  map1<-map(markers, tol = 10e-7)
  set_map_fun("haldane")
  list(map1,   marker_type(markers))
}

p <- c(2,2)
names(p) <- c("P1", "P2")
n.ind <- 1000
n.mrk <- 3
map.length <- 5
cm <- matrix(c("P1","P2"),
             ncol = 2, byrow = T)
D <- simulate_multiple_crosses(ploidy = p, #four parents
                               cross.mat = cm,
                               n.ind = rep(n.ind, nrow(cm)), # per cross
                               n.mrk= rep(n.mrk, length(p)), # per parent
                               alleles = list(c(1:4), c(1:4)),
                               map.length = map.length,
                               lambda = rep(2, length(p)))

D$phases
my_func(D)

(a1 <- D$map - min(D$map))
(a3 <- cumsum(mappoly::imf_h(c(0, my_func(D)[[1]]$seq.rf))))
H <- states_to_visit(D, err = 0.01, is.log = TRUE)
x1 <- est_map_R(states.hmm = H,tol = 10e-2)
(a2 <- cumsum(mappoly::imf_h(c(0, x1[[2]]))))
a2;a3


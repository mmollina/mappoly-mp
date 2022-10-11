require(mappoly2)
require(mappoly)
require(tidyverse)
require(foreach)
source("test_dev/sim_build_map.R")
ploidy.vec <- c(4, 2, 4, 6) #four parents
names(ploidy.vec) <- c("P1", "P2", "P3", "P4")
cross.mat <- matrix(c("P1","P2",
                      "P1","P3",
                      "P2","P3",
                      "P1","P1",
                      "P2","P4",
                      "P3","P4"), ncol = 2, byrow = T)
n.mrk <- c(100,100,100,100) #per parent
map.length <- 50 #in centimorgans

#### Parallel ####
n.cores <- parallel::detectCores()
cl <- parallel::makeCluster(n.cores)
parallel::clusterEvalQ(cl, require(mappoly2))
parallel::clusterEvalQ(cl, require(stringr))
parallel::clusterEvalQ(cl, require(mappoly))
doParallel::registerDoParallel(cl = cl)
n.sim <- 2000
#### completely informative, 100 ind per cross ####
system.time({
  alleles <- list(P1 = c(1:4),
                  P2 = c(5:6),
                  P3 = c(7:10),
                  P4 = c(11:12))
  n.ind <- c(100, 100, 100, 100, 100, 100) #per cross
  X1 <- foreach(i = 1:n.sim) %dopar% {
    myfunc(ploidy.vec,
           cross.mat,
           n.ind,
           n.mrk,
           alleles,
           map.length)
  }
}
)
#### biallelic, 100 ind per cross ####
system.time({
  alleles <- list(P1 = c(1:2),
                  P2 = c(1:2),
                  P3 = c(1:2),
                  P4 = c(1:2))
  n.ind <- c(100, 100, 100, 100, 100, 100) #per cross
  X2 <- foreach(i = 1:n.sim) %dopar% {
    myfunc(ploidy.vec,
           cross.mat,
           n.ind,
           n.mrk,
           alleles,
           map.length)
  }
}
)
#### completely informative, 30 ind per cross ####
system.time({
  alleles <- list(P1 = c(1:4),
                  P2 = c(5:6),
                  P3 = c(7:10),
                  P4 = c(11:12))
  n.ind <- c(30, 30, 30, 30, 30, 30) #per cross
  X3 <- foreach(i = 1:n.sim) %dopar% {
    myfunc(ploidy.vec,
           cross.mat,
           n.ind,
           n.mrk,
           alleles,
           map.length)
  }
}
)
#### biallelic, 30 ind per cross ####
system.time({
  alleles <- list(P1 = c(1:2),
                  P2 = c(1:2),
                  P3 = c(1:2),
                  P4 = c(1:2))
  n.ind <- c(30, 30, 30, 30, 30, 30) #per cross
  X4 <- foreach(i = 1:n.sim) %dopar% {
    myfunc(ploidy.vec,
           cross.mat,
           n.ind,
           n.mrk,
           alleles,
           map.length)
  }
}
)

#### Stop cluster ####
parallel::stopCluster(cl = cl)

#### Graphics ####
df1 <- reshape2::melt(X1, id.vars = c("map", "est.map", "mrks"), level = "rep")
df1 <- data.frame(df1, sim = "comp_info_100")
df2 <- reshape2::melt(X2, id.vars = c("map", "est.map", "mrks"), level = "rep")
df2 <- data.frame(df2, sim = "biallelic_100")
df3 <- reshape2::melt(X3, id.vars = c("map", "est.map", "mrks"), level = "rep")
df3 <- data.frame(df3, sim = "comp_info_30")
df4 <- reshape2::melt(X4, id.vars = c("map", "est.map", "mrks"), level = "rep")
df4 <- data.frame(df4, sim = "biallelic_30")
require(ggplot2)
DF <- rbind(df1, df2, df3, df4)
DF$Lrep <- as.factor(DF$Lrep)
save.image("test_dev/bias_study_multipop.rda")

ggplot(DF, aes(x = map, y = est.map) ) +
  stat_density_2d(aes(fill = stat(level)), geom = "polygon") +
  facet_grid(. ~ sim) + scale_fill_viridis_c() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")

ggplot(DF, aes(x = map, y = est.map, color = sim) ) +
  geom_point(alpha = .2) +
  geom_smooth(method = "lm", se = FALSE, lwd = 2) +
  facet_wrap(.~sim) + geom_abline(intercept = 0, slope = 1, linetype = "dashed")

ggplot(DF, aes(x = map, y = est.map) ) +
  geom_hex(bins = 50) + scale_fill_distiller(palette = "Spectral") +
  facet_wrap(.~sim) + geom_abline(intercept = 0, slope = 1, linetype = "dashed")



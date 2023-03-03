require(tidyverse)
g <- multi_fam_genoprob$genoprob

a <- unique(multi_fam_genoprob$pedigree)
a1 <- unlist(a[,1:2])
a2 <-  unlist(a[,3:4])
a <- unique(data.frame(a1, a2))

col <- NULL
n.found <- 10
pl <- c(4,4,4,4,4,4,4,4,4,4)
fc <- c('#e6194B',
        '#f58231',
        '#ffe119',
        '#3cb44b',
        '#469990',
        '#42d4f4',
        '#4363d8',
        '#000075',
        '#f032e6',
        '#000000')
cl<-fc[1:n.found]
for(i in 1:n.found) {
  fc <- colorRampPalette(c("white", cl[i], "black"))
  col <- c(col, fc(pl[i]*2)[(1+pl[i]/2):(pl[i]+pl[i]/2)])
}

barplot(rep(pl, each = 4), col = cl)
barplot(rep(pl, each = 4), col = col)


mp_pallet1 <- colorRampPalette(c("#ffe119", "#f58231","#e6194b","#808000","#9a6324", "#800000"))


Z1<- g %>% filter(ind == "X16400_N047")
head(Z1)
Z1 <- Z1 %>% mutate(parent = str_split_fixed(homolog, "_",2))

ggplot2::ggplot(Z1, ggplot2::aes(x = pos, y = prob, fill = homolog, color = homolog)) +
  ggplot2::geom_density(stat = "identity", alpha = 0.7, position = "stack") +
  ggplot2::facet_grid(rows = ggplot2::vars(homolog))

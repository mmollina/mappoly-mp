x<-est_map_R(H)
y<-choose(H$ploidy.cross.id[,1],H$ploidy.cross.id[,1]/2) * choose(H$ploidy.cross.id[,2],H$ploidy.cross.id[,2]/2)
w<-c(1, cumsum(y^2)+1)
i <- 1
z<-x[w[i]:(w[i+1]-1)]
dim(z) <- c(y[i],y[i])
image(z)
round(z,2)
dim(z)

# i: recombination fraction (rf[0]: between first and second)
# j: state in
# k: state out
# l: cross type

x<-est_map_R(H)
n_states<-choose(H$ploidy.cross.id[,1],H$ploidy.cross.id[,1]/2) * choose(H$ploidy.cross.id[,2],H$ploidy.cross.id[,2]/2)
accum_states<-c(1, cumsum(y^2)+1)

i <- 1

for(j in 1:4)
  for(k in 1:4)

x[i*n_states[l]*n_states[l] + j*n_states[l] + k + accum_states[l]*i]



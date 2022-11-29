M<-cor(dt[,n])
## eliminating NAs
id<-apply(M, 1, function(x) sum(is.na(x)) > length(n)*.8)
M <- abs(M[!id, !id])
dim(M)

heatmap(t(dt[,colnames(M[1:10, 1:10])]), Rowv = NA, Colv = NA)






round(M[1:10, 1:10], 2)
image(M)
A<-as.dist(1-M)
hc.snp <- hclust(A, method = "average")
dend.snp <- as.dendrogram(hc.snp)
plot(dend.snp)
groups.snp  <- cutree(tree = hc.snp, h = 0.05)
unique(groups.snp)




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




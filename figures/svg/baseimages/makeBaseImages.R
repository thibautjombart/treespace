library(treescape)
library(adegraphics)
library(ape)

## generate trees
set.seed(1)
all.D <- list()
for(i in 1:4) {
    all.D[[i]] <- dist(rnorm(10, rep(c(0,5), each=5)))
}
for(i in 5:7) {
    all.D[[i]] <- dist(1:10+runif(10))
}
for(i in 8:10) {
    all.D[[i]] <- dist(rep(1:2, c(5,5))+runif(10))
}

trees <- lapply(all.D, function(e) ladderize(root(nj(e),1)))
names(trees) <- paste("tree", 1:10, sep = "")
class(trees) <- "multiPhylo"
for(i in seq_along(trees)){
    trees[[i]]$tip.label <- letters[1:10]
}



## treescape result
res <- treescape(trees, method="patristic", nf=2)

## groves
groves <- findGroves(res, nclust=3)

## make graphics ##
## trees
for(i in seq_along(trees)){
    svg(paste("tree-",i,".svg", sep=""))
    plot(trees[[i]])
    dev.off()
}

## distances
svg("tablevalue.svg")
table.image(res$D, nclass=30)
dev.off()

## mds
svg("mds.svg")
plotGroves(res$pco, lab.show=TRUE, lab.cex=2)
dev.off()

## mds + clusters
svg("mdscolor.svg")
plotGroves(groves, lab.show=TRUE, col.pal=adegenet::seasun, lab.cex=2)
dev.off()

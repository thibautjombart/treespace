## ----setup, echo=FALSE---------------------------------------------------
# set global chunk options: images will be 7x5 inches
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.path="figs/", cache=FALSE)
options(digits = 4)

## ----load, message=FALSE-------------------------------------------------
library(treescape)

## ----tree1---------------------------------------------------------------
tree1 <- cbind(Infector=1:5,Infectee=2:6) 
tree1

## ----igraph_tree1, message=FALSE-----------------------------------------
library(igraph)
# set plotting options:
igraph_options(vertex.size=15,
               vertex.color="cyan",
               vertex.label.cex=2,
               edge.color="lightgrey",
               edge.arrow.size=1)

tree1graph <- graph_from_edgelist(tree1)
plot(tree1graph)

## ----simple_wiwMRCIs-----------------------------------------------------
findMRCIs(tree1)

## ----trees2_and_3--------------------------------------------------------
# a second scenario:
tree2 <- cbind(Infector=c(1,5,2,2,3),Infectee=2:6)
tree2
tree2graph <- graph_from_edgelist(tree2)
plot(tree2graph)

# and a third scenario:
tree3 <- cbind(Infector=c(2,2,2,2,6),Infectee=c(1,3,4,6,5)) 
tree3
tree3graph <- graph_from_edgelist(tree3)
plot(tree3graph)

## ----tree123_comparison--------------------------------------------------
m1 <- findMRCIs(tree1) # find the source case, MRCIs and MRCI depths for tree 1
m2 <- findMRCIs(tree2)
m3 <- findMRCIs(tree3)

matList <- list(m1$mrciDepths,m2$mrciDepths,m3$mrciDepths) # create a list of the mrciDepths matrices
matList
wiwTreeDist(matList, sampled=1:6) # find the Euclidean distances between these matrices, where all six cases are sampled

## ----tree123_sampled4:6--------------------------------------------------
wiwTreeDist(matList, sampled=4:6)

## ----trees1000-----------------------------------------------------------
set.seed(123)
num <- 500

# create a list of 500 random transmission trees with 11 cases, where the source case is fixed as case 1:
treelistSC1 <- lapply(1:num, function(x) {
  edges <- rtree(6)$edge # effectively creating a random transmission scenario
  relabel <- sample(1:11) # create a relabelling so that infections don't all happen in numerical order, but we force the source case to be 1:
  relabel[[which(relabel==1)]] <- relabel[[7]]
  relabel[[7]] <- 1
  relabelledEdges1 <- sapply(edges[,1], function(x) relabel[[x]])
  relabelledEdges2 <- sapply(edges[,2], function(x) relabel[[x]])
  cbind(relabelledEdges1,relabelledEdges2)
})

# create 500 more random transmission trees, but where the source case is fixed as case 2:
treelistSC2 <- lapply(1:num, function(x) {
  edges <- rtree(6)$edge 
  relabel <- sample(1:11) 
  relabel[[which(relabel==2)]] <- relabel[[7]]
  relabel[[7]] <- 2
  relabelledEdges1 <- sapply(edges[,1], function(x) relabel[[x]])
  relabelledEdges2 <- sapply(edges[,2], function(x) relabel[[x]])
  cbind(relabelledEdges1,relabelledEdges2)
})

# combine:
combinedLists <- c(treelistSC1,treelistSC2)

# get mrciDepths matrices:
matList1000 <- lapply(combinedLists, function(x)
  findMRCIs(x)$mrciDepths
)

# find pairwise tree distances, treating all cases as sampled:
WiwDists1000 <- wiwTreeDist(matList1000, sampled=1:11)

## ----wiw_MDS1000, message=FALSE------------------------------------------
wiwMDS <- dudi.pco(WiwDists1000, scannf=FALSE, nf=3)

library(ggplot2)
library(RColorBrewer)

wiwPlot <- ggplot(wiwMDS$li, aes(x=wiwMDS$li[,1],y=wiwMDS$li[,2]))

# prepare aesthetics
depths <- sapply(matList1000, function(x) mean(x))
sourcecase <- c(rep("1",num),rep("2",num))

# prepare colours:
colfunc <- colorRampPalette(brewer.pal(10,"Spectral"), space="Lab")

wiwPlot + 
  geom_point(size=4, colour="gray60", aes(shape=sourcecase)) + 
  geom_point(size=3, aes(colour=depths, shape=sourcecase)) +
  scale_colour_gradientn("Mean of v\n", 
                         colours=colfunc(7),
                         guide = guide_colourbar(barheight=10)) +
  scale_shape_discrete("Source case\n", solid=T, guide = guide_legend(keyheight = 3, keywidth=1.5)) +
  theme_bw(base_size = 12, base_family = "") +
  theme_bw(base_size = 12, base_family = "") +
  theme(
    legend.title = element_text(size=20),
    legend.text = element_text(size=20),
    axis.text.x = element_text(size=20), axis.text.y = element_text(size=20)) +
  xlab("") + ylab("")

## ----wiwMedian-----------------------------------------------------------
med <- wiwMedTree(matList1000)

## ----wiwMedian2----------------------------------------------------------
names(med)

## ----wiwMedTree----------------------------------------------------------
med$median

## ----wiwMedTreePlot------------------------------------------------------
medgraph <- graph_from_edgelist(combinedLists[[med$median]])
plot(medgraph)


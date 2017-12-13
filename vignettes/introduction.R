## ----setup, echo=FALSE---------------------------------------------------
# set global chunk options: images will be 7x5 inches
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.path="figs/", cache=FALSE, dpi=96)
options(digits = 4)
library("rgl")
knitr::knit_hooks$set(webgl=hook_webgl)

## ----install, eval=FALSE-------------------------------------------------
#  library(devtools)
#  install_github("thibautjombart/treespace")

## ----install2, eval=FALSE------------------------------------------------
#  install.packages("treespace")

## ----load----------------------------------------------------------------
library("treespace")

## ----load_packages, message=FALSE, warning=FALSE-------------------------
library("treespace")
library("adegenet")
library("adegraphics")
library("rgl")

## ----treespace-----------------------------------------------------------
# generate list of trees
set.seed(1)
x <- rmtree(10, 20)
names(x) <- paste("tree", 1:10, sep = "")

# use treespace
res <- treespace(x, nf=3)
names(res)
res

## ----distances-----------------------------------------------------------
# table.image
table.image(res$D, nclass=30)

# table.value with some customization
table.value(res$D, nclass=5, method="color", 
            symbol="circle", col=redpal(5))


## ----plotgroves----------------------------------------------------------
plotGroves(res$pco, lab.show=TRUE, lab.cex=1.5)

## ----plotgrovesD3--------------------------------------------------------
plotGrovesD3(res$pco, treeNames=1:10)

## ----woodmicePlots-------------------------------------------------------
data(woodmiceTrees)
wm.res <- treespace(woodmiceTrees,nf=3)

# PCs are stored in:
head(wm.res$pco$li)

# plot results
plotGrovesD3(wm.res$pco)

## ----findgroves, cache=FALSE---------------------------------------------
wm.groves <- findGroves(wm.res, nclust=6)
names(wm.groves)

## ----plotgroves2---------------------------------------------------------
# basic plot
plotGrovesD3(wm.groves)

# alternative with improved legend and tooltip text, giving the tree numbers:
plotGrovesD3(wm.groves, tooltip_text=paste0("Tree ",1:201), legend_width=50, col_lab="Cluster")

# plot axes 2 and 3. This helps to show why, for example, clusters 2 and 4 have been identified as separate, despite them appearing to overlap when viewing axes 1 and 2.
plotGrovesD3(wm.groves, xax=2, yax=3, tooltip_text=paste0("Tree ",1:201), legend_width=50, col_lab="Cluster")

## ----plotgroves_3D, rgl=TRUE, webgl=TRUE---------------------------------
# prepare a colour palette:
colours <- fac2col(wm.groves$groups, col.pal=funky)
plot3d(wm.groves$treespace$pco$li[,1],
       wm.groves$treespace$pco$li[,2],
       wm.groves$treespace$pco$li[,3],
       col=colours, type="s", size=1.5,
       xlab="", ylab="", zlab="")

## ----shiny_figures, echo=FALSE, fig.retina = NULL------------------------
knitr::include_graphics("figs/treespace3d.png", dpi=72)

knitr::include_graphics("figs/treespaceTree.png", dpi=72)

knitr::include_graphics("figs/treespaceDensiTree.png", dpi=72)

## ----woodmiceMedian------------------------------------------------------
# get first median tree
tre <- medTree(woodmiceTrees)$trees[[1]]

# plot tree
plot(tre,type="cladogram",edge.width=3, cex=0.8)

## ----woodmiceCluster1, out.width="600px"---------------------------------
# find median trees for the 6 clusters identified earlier:
res <- medTree(woodmiceTrees, wm.groves$groups)

# there is one output per cluster
names(res)

# get the first median of each
med.trees <- lapply(res, function(e) ladderize(e$trees[[1]]))

# plot trees
par(mfrow=c(2,3))
for(i in 1:length(med.trees)) plot(med.trees[[i]], main=paste("cluster",i),cex=1.5)


## ----woodmice_plotTreeDiff-----------------------------------------------
# Compare median trees from clusters 1 and 2:
plotTreeDiff(med.trees[[1]],med.trees[[2]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
# Compare median trees from clusters 1 and 4, and change aesthetics:
plotTreeDiff(med.trees[[1]],med.trees[[4]], type="cladogram", use.edge.length=FALSE, 
             treesFacing = TRUE, edge.width=2, colourMethod="palette", palette=spectral)

## ----woodmice-tip-emphasis-----------------------------------------------
wm3.res <- treespace(woodmiceTrees,nf=2,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)

# plot results
plotGrovesD3(wm3.res$pco)

## ----findgroves-with-emphasis--------------------------------------------
wm3.groves <- findGroves(woodmiceTrees,nf=3,nclust=6,emphasise.tips=c("No1007S","No1208S","No0909S"),emphasise.weight=3)
plotGrovesD3(wm3.groves)

## ----figure_construction, echo=FALSE, fig.retina = NULL------------------
knitr::include_graphics("figs/construction.png", dpi=72)

## ----treevec-------------------------------------------------------------
# generate a random tree:
tree <- rtree(6)
# topological vector of mrca distances from root:
treeVec(tree)
# vector of mrca distances from root when lambda=0.5:
treeVec(tree,0.5)
# vector of mrca distances as a function of lambda:
vecAsFunction <- treeVec(tree,return.lambda.function=TRUE)
# evaluate the vector at lambda=0.5:
vecAsFunction(0.5)

## ----treedist------------------------------------------------------------
# generate random trees
tree_a <- rtree(6)
tree_b <- rtree(6)

# topological (lambda=0) distance:
treeDist(tree_a,tree_b) 

# branch-length focused (lambda=1) distance:
treeDist(tree_a,tree_b,1)


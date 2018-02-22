## ----setup, echo=FALSE---------------------------------------------------
# set global chunk options: images will be 7x7 inches
knitr::opts_chunk$set(fig.width=7, fig.height=7, fig.path="figs/", cache=FALSE, dpi=96)
options(digits = 4)
library("rgl")
knitr::knit_hooks$set(webgl=hook_webgl)

## ----load, message=FALSE, warning=FALSE----------------------------------
library("treespace")
library("phangorn")
library("adegenet")
library("rgl")

## ----load_BEAST_trees----------------------------------------------------
data(DengueTrees)

## ----sample_BEAST_trees--------------------------------------------------
set.seed(123)
BEASTtrees <- DengueTrees[sample(1:length(DengueTrees),200)]

## ----load_seqs-----------------------------------------------------------
data(DengueSeqs)

## ----make_NJ-------------------------------------------------------------
makeTree <- function(x){
  tree <- nj(dist.dna(x, model = "TN93"))
  tree <- root(tree, resolve.root=TRUE, outgroup="D4Thai63")
  tree
}
DnjRooted <- makeTree(DengueSeqs)
# Note, there is a (small) negative branch length. 
# We set this to 0 to avoid warnings from the phangorn package later:
DnjRooted$edge.length[which(DnjRooted$edge.length < 0)] <- 0
plot(DnjRooted)

## ----make_NJ_boots, results="hide"---------------------------------------
Dnjboots <- boot.phylo(DnjRooted, DengueSeqs, B=100, 
	    	       makeTree, trees=TRUE, rooted=TRUE)
Dnjboots

## ----see_NJ_boots--------------------------------------------------------
plot(DnjRooted)
drawSupportOnEdges(Dnjboots$BP)

## ----make_ML, results="hide", message=FALSE------------------------------
Dfit.ini <- pml(DnjRooted, as.phyDat(DengueSeqs), k=4)
Dfit <- optim.pml(Dfit.ini, optNni=TRUE, optBf=TRUE,
                  optQ=TRUE, optGamma=TRUE, model="GTR")
# root:
DfitTreeRooted <- root(Dfit$tree, resolve.root=TRUE, outgroup="D4Thai63")

## ----view_ML-------------------------------------------------------------
plot(DfitTreeRooted)

## ----make_ML_boots, results="hide"---------------------------------------
# bootstrap supports:
DMLboots <- bootstrap.pml(Dfit, optNni=TRUE)
# root:
DMLbootsrooted <- lapply(DMLboots, function(x) root(x, resolve.root=TRUE, outgroup="D4Thai63"))
class(DMLbootsrooted) <- "multiPhylo"


## ----see_ML_boots--------------------------------------------------------
plotBS(DfitTreeRooted, DMLboots, type="phylogram")

## ----run_treespace-------------------------------------------------------
# collect the trees into a single object of class multiPhylo:
DengueTrees <- c(BEASTtrees, Dnjboots$trees, DMLbootsrooted,
		             DnjRooted, DfitTreeRooted)
class(DengueTrees) <- "multiPhylo"
# add tree names:
names(DengueTrees)[1:200] <- paste0("BEAST",1:200)
names(DengueTrees)[201:300] <- paste0("NJ_boots",1:100)
names(DengueTrees)[301:400] <- paste0("ML_boots",1:100)
names(DengueTrees)[[401]] <- "NJ"
names(DengueTrees)[[402]] <- "ML"
# create vector corresponding to tree inference method:
Dtype <- c(rep("BEAST",200),rep("NJboots",100),rep("MLboots",100),"NJ","ML")

# use treespace to find and project the distances:
Dscape <- treespace(DengueTrees, nf=5)

## ----simple_plot---------------------------------------------------------
# simple plot:
plotGrovesD3(Dscape$pco, groups=Dtype)

## ----make_better_plot----------------------------------------------------
Dcols <- c("#1b9e77","#d95f02","#7570b3")
Dmethod <- c(rep("BEAST",200),rep("NJ",100),rep("ML",100),"NJ","ML")
Dbootstraps <- c(rep("replicates",400),"NJ","ML")
Dhighlight <- c(rep(1,400),2,2)
plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             colors=Dcols,
             col_lab="Tree type",
             size_var=Dhighlight,
             size_range = c(100,500),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,400),1,1), 
             legend_width=80)

## ----make_better_plot_with_labels----------------------------------------
plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             treeNames = names(DengueTrees), # add the tree names as labels
             colors=Dcols,
             col_lab="Tree type",
             size_var=Dhighlight,
             size_range = c(100,500),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,400),1,1), 
             legend_width=80)

## ----make_better_plot_with_tooltips--------------------------------------
plotGrovesD3(Dscape$pco, 
             groups=Dmethod, 
             tooltip_text = names(DengueTrees), # add the tree names as tooltip text
             colors=Dcols,
             col_lab="Tree type",
             size_var=Dhighlight,
             size_range = c(100,500),
             size_lab="",
             symbol_var=Dbootstraps,
             symbol_lab="",
             point_opacity=c(rep(0.4,400),1,1), 
             legend_width=80)

## ----scree_plot----------------------------------------------------------
barplot(Dscape$pco$eig, col="navy")

## ----load_rgl------------------------------------------------------------
library(rgl)

## ----plot_3D, rgl=TRUE, webgl=TRUE---------------------------------------
Dcols3D <- c(rep(Dcols[[1]],200),rep(Dcols[[2]],100),rep(Dcols[[3]],100),Dcols[[2]],Dcols[[3]])
rgl::plot3d(Dscape$pco$li[,1],Dscape$pco$li[,2],Dscape$pco$li[,3],
       type="s",
       size=c(rep(1.5,400),3,3), 
       col=Dcols3D,
       xlab="", ylab="", zlab="")

## ----NJ_and_ML_overlap---------------------------------------------------
# trees with the same topology as the NJ tree:
which(as.matrix(Dscape$D)["NJ",]==0)
# trees with the same topology as the ML tree:
which(as.matrix(Dscape$D)["ML",]==0)

## ----compare_trees_NJ_v_ML-----------------------------------------------
# comparing NJ and ML:
plotTreeDiff(DnjRooted,DfitTreeRooted, use.edge.length=FALSE)
treeDist(DnjRooted,DfitTreeRooted)

## ----compare_trees_NJ_v_ML_recoloured------------------------------------
# comparing NJ and ML:
plotTreeDiff(DnjRooted,DfitTreeRooted, use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)

## ----make_BEAST_median---------------------------------------------------
BEASTmed <- medTree(BEASTtrees)

## ----compare_BEAST_meds--------------------------------------------------
BEASTmed$trees
treeDist(BEASTmed$trees[[1]],BEASTmed$trees[[2]])

## ----save_BEAST_median---------------------------------------------------
BEASTrep <- BEASTmed$trees[[1]]

## ----compare_BEAST_to_other_trees----------------------------------------
# comparing BEAST median and NJ:
plotTreeDiff(BEASTrep,DnjRooted, use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
treeDist(BEASTrep,DnjRooted)
# comparing BEAST median and ML:
plotTreeDiff(BEASTrep,DfitTreeRooted, use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
treeDist(BEASTrep,DfitTreeRooted)
# comparing BEAST median to a random BEAST tree:
num <- runif(1,1,200)
randomBEASTtree <- BEASTtrees[[num]]
plotTreeDiff(BEASTrep, randomBEASTtree, use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
treeDist(BEASTrep,randomBEASTtree)

## ----BEASTtrees----------------------------------------------------------
# load the MCC tree
data(DengueBEASTMCC)
# concatenate with other BEAST trees
BEAST201 <- c(BEASTtrees, DengueBEASTMCC)
# compare using treespace:
BEASTscape <- treespace(BEAST201, nf=5)
# simple plot:
plotGrovesD3(BEASTscape$pco)

## ----BEASTtrees_clusters-------------------------------------------------
# find clusters or 'groves':
BEASTGroves <- findGroves(BEASTscape, nclust=4, clustering = "single")

## ----BEASTtrees_meds-----------------------------------------------------
# find median tree(s) per cluster:
BEASTMeds <- medTree(BEAST201, groups=BEASTGroves$groups)
# for each cluster, select a single median tree to represent it:
BEASTMedTrees <- c(BEASTMeds$`1`$trees[[1]],
                   BEASTMeds$`2`$trees[[1]],
                   BEASTMeds$`3`$trees[[1]],
                   BEASTMeds$`4`$trees[[1]])

## ----BEASTtrees_plot, warning=FALSE--------------------------------------
# extract the numbers from the tree list 'BEASTtrees' which correspond to the median trees: 
BEASTMedTreeNums <-c(which(BEASTGroves$groups==1)[[BEASTMeds$`1`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==2)[[BEASTMeds$`2`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==3)[[BEASTMeds$`3`$treenumbers[[1]]]],
                     which(BEASTGroves$groups==4)[[BEASTMeds$`4`$treenumbers[[1]]]])
# prepare a vector to highlight median and MCC trees
highlightTrees <- rep(1,201)
highlightTrees[[201]] <- 2
highlightTrees[BEASTMedTreeNums] <- 2
# prepare colours:
BEASTcols <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3")

# plot:
plotGrovesD3(BEASTscape$pco,
          groups=as.vector(BEASTGroves$groups),
          colors=BEASTcols,
          col_lab="Cluster",
          symbol_var = highlightTrees,
          size_range = c(60,600),
          size_var = highlightTrees,
          legend_width=0)

## ----BEASTtree_diffs-----------------------------------------------------
# differences between the MCC tree and the median from the largest cluster:
treeDist(DengueBEASTMCC,BEASTMedTrees[[1]])
plotTreeDiff(DengueBEASTMCC,BEASTMedTrees[[1]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)
# differences between the median trees from clusters 1 and 2:
treeDist(BEASTMedTrees[[1]],BEASTMedTrees[[2]])
plotTreeDiff(BEASTMedTrees[[1]],BEASTMedTrees[[2]], use.edge.length=FALSE, 
             treesFacing = TRUE, colourMethod = "palette", palette = funky)


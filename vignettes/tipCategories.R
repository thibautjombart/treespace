## ----setup, echo=FALSE--------------------------------------------------------
# set global chunk options: images will be 7x4 inches
knitr::opts_chunk$set(fig.width=7, fig.height=4, cache=FALSE, dpi=96)
options(digits = 4)

## ----setupVisible, message=FALSE----------------------------------------------
# load treespace and packages for plotting:
library(treespace)
library(RColorBrewer) 
library(ggplot2)
library(reshape2)
# set colour scheme
pal <- brewer.pal(3,"Dark2")

## ----create_examples_table, echo=FALSE----------------------------------------
cats <- c("Bacterial sub-types e.g. serogroups",
          "Species",
          "Host",
          "Protein families",
          "Population groups",
          "Disjoint features or phenotypes",
          "Broad taxonomy")
indivs <- c("Bacterial isolates",
            "Orthologous genes",
            "Deep sequencing reads of pathogen",
            "Proteins",
            "Individual organisms",
            "Individual organisms",
            "Individual organisms")
exTable <- cbind(cats,indivs)
colnames(exTable) <- c("Categories", "Individuals")

## ----table1, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'--------
require(pander)
panderOptions('table.split.table', Inf)
pander(exTable, style = 'rmarkdown')

## ----create_trees_for_collapsing_example, echo=FALSE--------------------------
tr1 <- read.tree(text="(((c4,c3),(c2,c1)),((b3,(b2,b1)),((a3,a2),a1)));")
tr1$tip.label <- c("Patient C read 4","Patient C read 3",
                  "Patient C read 2","Patient C read 1",
                  "Patient B read 3","Patient B read 2",
                  "Patient B read 1",
                  "Patient A read 3","Patient A read 2",
                  "Patient A read 1")

tr1Collapsed <- read.tree(text="(C,(B,A));")
tr1Collapsed$tip.label <- c("Patient C","Patient B","Patient A")


tr2 <- read.tree(text="((((c4,c3),(c2,c1)),b2),(b1,((a3,a2),a1)));")
tr2$tip.label <- c("Patient C read 4","Patient C read 3",
                  "Patient C read 2","Patient C read 1",
                  "Patient B read 2","Patient B read 1",
                  "Patient A read 3","Patient A read 2",
                  "Patient A read 1")

tr2Collapsed <- read.tree(text="((C,B),(B,A));")
tr2Collapsed$tip.label <- c("Patient C","Patient B","Patient B","Patient A")

## ----plot_patient_trees, echo=FALSE-------------------------------------------
layout(matrix(1:2,1,2))
plot(tr1, tip.color=c(rep(pal[[1]],4),rep(pal[[2]],3),rep(pal[[3]],3)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, cex=0.8,
     edge.color=c("black",rep(pal[[1]],6),"black",rep(pal[[2]],5),rep(pal[[3]],5)))
plot(tr2, tip.color=c(rep(pal[[1]],4),rep(pal[[2]],2),rep(pal[[3]],3)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, cex=0.8,
     edge.color=c(rep("black",2),rep(pal[[1]],6),pal[[2]],"black",pal[[2]],"black",      rep(pal[[3]],4)))

## ----plot_collapsed_trees, echo=FALSE-----------------------------------------
layout(matrix(1:2,1,2))
plot(tr1Collapsed, tip.color=c(pal[[1]],pal[[2]],pal[[3]]),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset=0.1, font=4, cex=0.8,
     edge.color=c(pal[[1]],"black",pal[[2]],pal[[3]]))
plot(tr2Collapsed, tip.color=c(pal[[1]],rep(pal[[2]],2),pal[[3]]),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset=0.1, font=4, cex=0.8,
     edge.color=c("black",pal[[1]],pal[[2]],"black",pal[[2]],pal[[3]]))

## ----relatedTreeDist----------------------------------------------------------
df <- cbind(c(rep("Patient A",3),rep("Patient B",3),rep("Patient C",4)),
            sort(tr1$tip.label))
df
relatedTreeDist(list(tr1,tr2),df)[[1]]

## ----tipsMRCAdepths-----------------------------------------------------------
tipsMRCAdepths(tr1Collapsed)
tipsMRCAdepths(tr2Collapsed)

## ----calculation_of_relatedTreeDist, echo=FALSE-------------------------------
MRCAdepths <- cbind(tipsMRCAdepths(tr1Collapsed),c(0.5,0,0.5))
colnames(MRCAdepths) <- c("tip1","tip2","Tree 1", "Tree 2")
MRCAdepths

## ----six_comparable_trees, fig.height=12--------------------------------------
suppressWarnings(RNGversion("3.5.0"))
set.seed(948)
# set colour scheme
pal2 <- brewer.pal(8,"Dark2")

# create a "base" (category-level) tree
baseTree <- rtree(8)
baseTree$tip.label <- letters[8:1]

tree1 <- simulateIndTree(baseTree, itips=3, permuteTips=FALSE)
tree2 <- simulateIndTree(baseTree, itips=4, permuteTips=FALSE)
tree2$tip.label <- c(paste0("h_",2:5),paste0("g_",3:6),paste0("f_",c(2,3,7,9)),
                     paste0("e_",c(1,2,5,6)),paste0("d_",3:6),paste0("c_",5:8),
                     paste0("b_",c(3,5,6,9)),paste0("a_",c(1,4,8,9)))
tree3 <- simulateIndTree(baseTree, itips=4, tipPercent = 20)
tree3NotForPlotting <- simulateIndTree(baseTree, itips=4, permuteTips=FALSE) # just for setting colours later 
tree4 <- simulateIndTree(baseTree, itips=6)
tree4NotForPlotting <- simulateIndTree(baseTree, itips=6, permuteTips=FALSE)

# create another base tree
baseTree2 <- rtree(8, tip.label=letters[8:1])
tree5 <- simulateIndTree(baseTree2, itips=6, permuteTips=FALSE)
tree6 <- simulateIndTree(baseTree2, itips=6, tipPercent=30)

# set up colour palettes
tipcolors3 <- c(rep(pal2[[1]],3),rep(pal2[[2]],3),rep(pal2[[3]],3),rep(pal2[[4]],3),rep(pal2[[5]],3),rep(pal2[[6]],3),rep(pal2[[7]],3),rep(pal2[[8]],3))
tipcolors4 <- c(rep(pal2[[1]],4),rep(pal2[[2]],4),rep(pal2[[3]],4),rep(pal2[[4]],4),rep(pal2[[5]],4),rep(pal2[[6]],4),rep(pal2[[7]],4),rep(pal2[[8]],4)) # colours for 4 tips
tipcolors6 <- c(rep(pal2[[1]],6),rep(pal2[[2]],6),rep(pal2[[3]],6),rep(pal2[[4]],6),rep(pal2[[5]],6),rep(pal2[[6]],6),rep(pal2[[7]],6),rep(pal2[[8]],6)) # colours for 6 tips

# prepare tip colours for plotting
tree3TipOrder <- sapply(tree3$tip.label, function(x) which(tree3NotForPlotting$tip.label==x))
tree4TipOrder <- sapply(tree4$tip.label, function(x) which(tree4NotForPlotting$tip.label==x))
tree5TipOrder <- sapply(tree5$tip.label, function(x) which(tree4NotForPlotting$tip.label==x))
tree6TipOrder <- sapply(tree6$tip.label, function(x) which(tree4NotForPlotting$tip.label==x))

layout(matrix(c(1,4,2,5,3,6), 2,3))
plot(tree1, tip.color=tipcolors3, no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=2)
plot(tree2, tip.color=tipcolors4, no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=2)
plot(tree3, tip.color=tipcolors4[tree3TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.8)
plot(tree4, tip.color=tipcolors6[tree4TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)
plot(tree5, tip.color=tipcolors6[tree5TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)
plot(tree6, tip.color=tipcolors6[tree6TipOrder], no.margin=TRUE,
     edge.width = 4, use.edge.length = FALSE,
     label.offset= 0.5, font=4, cex=1.2)

## ----relatedTreeDist_six_trees------------------------------------------------
trees <- list(tree1,tree2,tree3,tree4,tree5,tree6)
df <- cbind(sort(rep(letters[1:8],9)),sort(paste0(letters[1:8],"_",rep(1:9,8))))

dists <- relatedTreeDist(trees,df)
dists

## ----six_heatmap--------------------------------------------------------------
dists <- as.matrix(dists)
colnames(dists) <- rownames(dists) <- c("Tree 1", "Tree 2", "Tree 3", "Tree 4",
                                        "Tree 5", "Tree 6")

melted_dists <- melt(dists, na.rm=TRUE)

ggheatmap <- ggplot(data = melted_dists, aes(Var2, Var1, fill = value))+
  geom_tile(color = "darkgrey")+
  scale_fill_gradient2(low = "white", high = "firebrick2",
                       name="Tree distance") +
  theme_minimal() + coord_fixed()

ggheatmap +
  geom_text(aes(Var2, Var1, label = signif(value,2)), color = "black", size = 8) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text = element_text(size=12),
    axis.ticks = element_blank(),
    legend.position = "none" )


## ----concordance_basic--------------------------------------------------------
catTree <- read.tree(text="(C,(B,A));")
indTree1 <- read.tree(text="(((c4,c3),(c2,c1)),((b1,b2),((a3,a2),a1)));")
indTree2 <- read.tree(text="(((c4,c3),(c2,c1)),((b1,a2),((a3,b2),a1)));")
indTree3 <- read.tree(text="((a3,(a2,a1)),((b1,c2),((c3,b2),(c1,c4))));")

plot(catTree, tip.color=pal,
     edge.width = 4, type="cladogram",
     label.offset= 0.5, font=4,
     edge.color=c(pal[[1]],"black",pal[[2]],pal[[3]]))
layout(matrix(1:3,1,3))
plot(indTree1, tip.color=c(rep(pal[[1]],4),rep(pal[[2]],2),rep(pal[[3]],3)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[1]],6),"black",rep(pal[[2]],3),rep(pal[[3]],5)))
plot(indTree2, tip.color=c(rep(pal[[1]],4),pal[[2]],rep(pal[[3]],2),pal[[2]],pal[[3]]),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[1]],6),rep("black",2),pal[[2]],pal[[3]],
                  rep("black",2),pal[[3]],pal[[2]],pal[[3]]))
plot(indTree3, tip.color=c(rep(pal[[3]],3),pal[[2]],rep(pal[[1]],2),pal[[2]],rep(pal[[1]],2)),
     edge.width = 4, type="cladogram", no.margin=TRUE,
     label.offset= 0.5, font=4, cex=2,
     edge.color=c("black",rep(pal[[3]],4),rep("black",2),pal[[2]],pal[[1]],
                  rep("black",2),pal[[1]],pal[[2]],rep(pal[[1]],3)))


## ----concordance_T1-----------------------------------------------------------
df <- cbind(c(rep("A",3),rep("B",2),rep("C",4)),sort(indTree1$tip.label))
treeConcordance(catTree,indTree1,df)

## ----concordance_T2-----------------------------------------------------------
treeConcordance(catTree,indTree2,df)

## ----concordance_T3-----------------------------------------------------------
treeConcordance(catTree,indTree3,df)

## ----concordance_permuations--------------------------------------------------
n <- 5
reps <- 10
reftree <- rtree(n, tip.label=letters[1:n])
indTrees <- lapply(rep(seq(0,100,20),reps), function(x)
  simulateIndTree(reftree,itips=n,permuteTips=TRUE,tipPercent=x))

df <- cbind(sort(rep(letters[1:n],n)),sort(indTrees[[1]]$tip.label))

concordances <- sapply(indTrees, function(x) treeConcordance(reftree,x,df))

resultsTab <- as.data.frame(cbind(rep(seq(0,100,20),reps),concordances))
colnames(resultsTab) <- c("Percentage","Concordance")
resultsTab[,1] <- factor(resultsTab[,1], levels=seq(0,100,20))
plot <- ggplot(resultsTab, aes(x=Percentage, y=Concordance))

plot + geom_boxplot(aes(colour=Percentage)) + theme_bw() + guides(colour=FALSE) +
  xlab("Percentage of tips permuted") + ylim(c(0,1)) +
  theme(axis.text = element_text(size=18),
        axis.title = element_text(size=18))


#' Find tip position differences
#'
#' Find the topologicial differences between two trees with the same tip labels. The function returns a data frame of the tips and the number of differences in their ancestry between the two trees.
#' Called by \code{\link{plotTreeDiff}}, which highlights the differing tips in a plot of the two trees.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom combinat combn2
#'
#' @param tr1 an object of the class \code{phylo}: the first tree to compare.
#' @param tr2 an object of the class \code{phylo}: the second tree to compare.
#' @param vec1 an optional input, the result of \code{treeVec(tr1, lambda=0)}, to speed up the computation.
#' @param vec2 an optional input, the result of \code{treeVec(tr2, lambda=0)}, to speed up the computation.
#'
#' @return
#' A data frame of the tree tips and the number of ancestral differences between them in the two trees, in order of increasing difference.
#' A tip is said to have zero difference if each of its ancestral nodes admits the same tip partition in the two trees.
#' 
#' @seealso \code{\link{medTree}} \code{\link{plotTreeDiff}}
#'
#' @examples
#' ## simple example on trees with five tips:
#' tr1 <- read.tree(text="((A:1,B:1):1,((C:1,D:1):1,E:1):1):1;")
#' tr2 <- read.tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1):1;")
#' tipDiff(tr1,tr2)
#' 
#' ## example on larger woodmice trees
#' data(woodmiceTrees)
#' tipDiff(woodmiceTrees[[1]],woodmiceTrees[[2]])
#' 
#' @export 
tipDiff <- function(tr1,tr2,vec1=NULL,vec2=NULL) {
  
  l <- length(tr1$tip.label)
  lchoose2 <- l*(l-1)/2
  if( l != length(tr2$tip.label)) stop("Trees must have the same number of tips")
  
  if(setequal(tr1$tip.label,tr2$tip.label) == FALSE) stop("Trees must have the same tip label sets")
  
  # get vec1
  if (is.null(vec1)) {
    vec1 <- treeVec(tr1) # emphasise.tips, emphasise.weight)
  }
  if (is.null(vec2)) {
    vec2 <- treeVec(tr2) # emphasise.tips, emphasise.weight)
  }
  
  # trim pendant edge entries from vectors
  vec1 <- vec1[1:lchoose2] # emphasise.tips, emphasise.weight)
  vec2 <- vec2[1:lchoose2]
  
  # find the positions where the vectors are different
  vecDiff <- (vec1!=vec2)
  
  # combine the tip labels (in order) and whether the vectors are different
  treedf <- as.data.frame(cbind(combn2(tr1$tip.label[order(tr1$tip.label)]),vecDiff))
  
  # list the tips which appear in "TRUE" rows
  allTipDiffs <- c(as.character(treedf[,1][vecDiff]),as.character(treedf[,2][vecDiff]))
  
  # find the number of times each tip appears
  tipSignificance <- sapply(tr1$tip.label, function(x)
    length(which(allTipDiffs==x)))
  
  # prepare output as data frame of tips and their differences
  tipDiff <- as.data.frame(cbind(names(tipSignificance[order(tipSignificance)]),tipSignificance[order(tipSignificance)]), stringsAsFactors = FALSE)
  class(tipDiff[,2]) <- "numeric"
  rownames(tipDiff) <- NULL
  colnames(tipDiff) <- c("Tip","No. of differences")
  
  return(tipDiff)
}




#' Plot tree differences
#'
#' Highlight the topologicial differences between two trees, plotted side by side. 
#' This function is useful for comparing representative "median" trees - see \code{\link{medTree}}.
#' It relies on the function \code{\link{tipDiff}}.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom adegenet lightseasun
#' @importFrom adegenet num2col
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics layout
#'
#' @param tr1 an object of the class \code{phylo}: the first tree to plot.
#' @param tr2 an object of the class \code{phylo}: the second tree to plot.
#' @param tipDiff an optional input, the result of \code{\link{tipDiff}}. Supplying this will save time if calling \code{plotTreeDiff} repeatedly, for example with different aesthetics.
#' @param vec1 an optional input, the result of \code{treeVec(tr1, lambda=0)}. This argument is ignored if \code{tipDiff} is supplied; otherwise supplying this will save time if calling \code{plotTreeDiff} repeatedly, for example with different aesthetics.
#' @param vec2 an optional input, the result of \code{treeVec(tr2, lambda=0)}. This argument is ignored if \code{tipDiff} is supplied; otherwise supplying this will save time if calling \code{plotTreeDiff} repeatedly, for example with different aesthetics.
#' @param baseCol the colour used for tips with identical ancestry in the two trees. Defaults to "grey".
#' @param colourMethod the method to use for colouring. Default is "ramp", corresponding to the original implementation, where the function \code{colorRampPalette} is used to create a palette which ranges from \code{col1} to \code{col2}. For large trees this can be hard to interpret, and method \code{palette} may be preferred, which permits the selection of a palette to use in \code{adegenet}'s function \code{num2col}.
#' @param col1 the first colour used to define the colour spectrum for tips with differences. This colour will be used for tips with minor differences. Defaults to "peachpuff". Ignored if \code{colourMethod="palette"}.
#' @param col2 the second colour used to define the colour spectrum for tips with differences. This colour will be used for tips with major differences. Defaults to "red2". Ignored if \code{colourMethod="palette"}.
#' @param palette the colour palette to be used if \code{colourMethod="palette"}. For a list of available palettes see \code{?num2col}.
#' @param ... further arguments passed to \code{\link{plot.phylo}}
#'
#' @return
#' A plot of the two trees side by side. Tips are coloured in the following way:
#' \itemize{
#' \item if each ancestor of a tip in tree 1 occurs in tree 2 with the same partition of tip descendants, then the tip is coloured grey (or supplied "baseCol")
#' \item if not, the tip gets coloured pale orange to red on a scale according to how many differences there are amongst its most recent common ancestors with other tips. The colour spectrum can be changed according to preference.
#' }
#' 
#' @seealso \code{\link{medTree}}, \code{\link{tipDiff}}
#'
#' @examples
#' ## simple example on trees with five tips:
#' tr1 <- read.tree(text="((A:1,B:1):1,((C:1,D:1):1,E:1):1):1;")
#' tr2 <- read.tree(text="((A:1,B:1):1,(C:1,(D:1,E:1):1):1):1;")
#' plotTreeDiff(tr1,tr2)
#' 
#' ## example on larger woodmice trees
#' data(woodmiceTrees)
#' # find the tip differences in advance, to avoid recalculating with each plot
#' wmTipDiff <- tipDiff(woodmiceTrees[[1]],woodmiceTrees[[2]])
#' plotTreeDiff(woodmiceTrees[[1]],woodmiceTrees[[2]], tipDiff=wmTipDiff)
#' ## change aesthetics:
#' plotTreeDiff(woodmiceTrees[[1]],woodmiceTrees[[2]], tipDiff=wmTipDiff,
#'    baseCol="grey2", col1="cyan", col2="navy", 
#'    edge.width=2, type="radial", cex=0.5, font=2)
#' ## use colour palette from adegenet:
#' plotTreeDiff(woodmiceTrees[[1]],woodmiceTrees[[2]], tipDiff=wmTipDiff,
#'    baseCol="black", colourMethod="palette", 
#'    edge.width=2, type="cladogram", cex=0.5, font=2)    
#' 
#' @export 
plotTreeDiff <- function(tr1,tr2,tipDiff=NULL,vec1=NULL,vec2=NULL,baseCol="grey",col1="peachpuff",col2="red2",colourMethod="ramp",palette=lightseasun,...) {
  
  l <- length(tr1$tip.label)
  
  if (is.null(tipDiff)){
    # call tipDiff
    tipDiff <- tipDiff(tr1,tr2,vec1,vec2)
  }
  
  # find the number of times each tip appears, in the order that the tip labels are read
  tipSignificance1 <- sapply(tr1$tip.label, function(x)
    tipDiff[which(tipDiff[,1]==x),2])
  tipSignificance2 <- sapply(tr2$tip.label, function(x)
    tipDiff[which(tipDiff[,1]==x),2])
  
  if (colourMethod=="ramp") {
    colfunc <- colorRampPalette(c(col1,col2))
  
    if (min(tipDiff[,2])==0) { # make sure tips with no differences are coloured baseCol
      tipSignificance1 <- tipSignificance1 + 1
      tipSignificance2 <- tipSignificance2 + 1
      if (max(tipDiff[,2])==0) {numCols <- 0}
      else {numCols <- max(tipDiff[,2]) - min(tipDiff[,2][which(tipDiff[,2]!=0)]) + 1}
      pal <- c(baseCol,colfunc(numCols))
    }
    else {
      numCols <- max(tipDiff[,2]) - min(tipDiff[,2]) + 1
      pal <- colfunc(numCols)
    }
  
    tipCols1 <- pal[as.factor(tipSignificance1)]
    tipCols2 <- pal[as.factor(tipSignificance2)]
  }
  
  else {
    tipCols1 <- num2col(tipSignificance1, col.pal=palette)
    tipCols2 <- num2col(tipSignificance2, col.pal=palette)
      if (min(tipDiff[,2])==0) { # make sure tips with no differences are coloured baseCol
        tipCols1[which(tipSignificance1==0)] <- baseCol
        tipCols2[which(tipSignificance2==0)] <- baseCol
      }
  }
    
  # plot
  layout(matrix(1:2, 1, 2))
  plot.phylo(tr1, tip.color=tipCols1, no.margin=TRUE, ...)
  plot.phylo(tr2, tip.color=tipCols2, no.margin=TRUE, ...)
  
}
#' Find tip position differences
#'
#' Find the topologicial differences between two trees with the same tip labels. The function returns a data frame of the tips and the number of differences in their ancestry between the two trees.
#' Called by \code{\link{plotTreeDiff}}, which highlights the differing tips in a plot of the two trees.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom combinat combn2
#' @importFrom adegenet num2col
#'
#' @param tr1 an object of the class \code{phylo}: the first tree to compare.
#' @param tr2 an object of the class \code{phylo}: the second tree to compare.
#' @param vec1 an optional input, the result of \code{treeVec(tr1, lambda=0)}, to speed up the computation.
#' @param vec2 an optional input, the result of \code{treeVec(tr2, lambda=0)}, to speed up the computation.
#' @param sizeOfDifferences a logical (default FALSE) specifying whether the size of the differences in the vectors per tip is also computed
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
tipDiff <- function(tr1,tr2,vec1=NULL,vec2=NULL,sizeOfDifferences=FALSE) {
  
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
  orderedLabels <- tr1$tip.label[order(tr1$tip.label)]
  treedf <- cbind.data.frame(combn2(orderedLabels),vecDiff)
  
  # extract the "TRUE" rows
  treeDiffsDf <- treedf[vecDiff,]
  tipsWithDifferencesRepeated <- unlist(list(treeDiffsDf[,1],treeDiffsDf[,2])) # list the tips with differences, possibly (probably) with repeats
  tipsWithDifferencesNames <- unique(tipsWithDifferencesRepeated) # get the names of the tips with differences (each listed just once)
  
  # find the number of times each tip appears
  tipSignificance <- sapply(orderedLabels, function(x)
    if (x %in% tipsWithDifferencesNames) length(which(tipsWithDifferencesRepeated==x))
    else 0
  )
  
  # prepare output as data frame of tips and their differences
  orderOfSignificance <- order(tipSignificance)
  tipDiff <- cbind.data.frame(orderedLabels[orderOfSignificance],tipSignificance[orderOfSignificance])
  rownames(tipDiff) <- NULL
  colnames(tipDiff) <- c("Tip","No. of differences")
  
  
  if (sizeOfDifferences) { # also compute the sizes of the differences
    treedfSizes <- cbind.data.frame(treedf[vecDiff,],abs(vec1[vecDiff]-vec2[vecDiff]))
    
    tipSizeOfDifferences <- sapply(orderedLabels[orderOfSignificance], function(x)
      sum(treedfSizes[which(treedfSizes[,1]==x),4]) + sum(treedfSizes[which(treedfSizes[,2]==x),4])
    )
    
    tipDiff <- cbind.data.frame(tipDiff,tipSizeOfDifferences)
    rownames(tipDiff) <- NULL
    colnames(tipDiff) <- c("Tip","No. of differences","Size of differences")
  }
    
  return(tipDiff)
}




#' Plot tree differences
#'
#' Highlight the topological differences between two trees, plotted side by side. 
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
#' @param sizeOfDifferences a logical (default FALSE) specifying whether the size of the tip differences should be used, or just a count of the number of differences (see \code{tipDiff})
#' @param tipMatch a logical (default TRUE) specifying whether the second tree should be rotated so that, as far as possible, each of its tips lies opposite its equivalent in the first tree
#' @param treesFacing a logical (default FALSE) specifying whether the trees should be plotted facing each other - that is, with the second tree plotted "leftwards".
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
#' tr1 <- woodmiceTrees[[1]]
#' tr2 <- woodmiceTrees[[57]] # for example
#' 
#' # find the tip differences in advance, to avoid recalculating with each plot
#' wmTipDiff <- tipDiff(tr1,tr2, sizeOfDifferences=TRUE)
#' plotTreeDiff(tr1,tr2, tipDiff=wmTipDiff, tipMatch=TRUE)
#' 
#' ## change aesthetics:
#' # trees facing each other:
#' plotTreeDiff(tr1,tr2, tipDiff=wmTipDiff, treesFacing=TRUE)
#' 
#' # radial plots, and change colours:
#' plotTreeDiff(tr1,tr2, tipDiff=wmTipDiff,
#'    baseCol="grey2", col1="cyan", col2="navy", 
#'    edge.width=2, type="radial", cex=0.5, font=2)
#' # cladogram plots, and use colour palette from adegenet to see differences more clearly:
#' plotTreeDiff(tr1,tr2, tipDiff=wmTipDiff,
#'    treesFacing=TRUE, baseCol="black", colourMethod="palette", 
#'    edge.width=2, type="cladogram", cex=0.5, font=2)    
#'
#' # including the size of the differences highlights tip "No0906s" a little more:
#' # (this is typically a more informative plot in cases where many tips have the 
#' # same difference count, for example when a whole clade has been shifted "up" 
#' # or "down" the tree but its internal topology remains the same.) 
#'  
#' plotTreeDiff(tr1,tr2, tipDiff=wmTipDiff, sizeOfDifferences=TRUE,
#'    treesFacing=TRUE, baseCol="black", colourMethod="palette", 
#'    edge.width=2, type="cladogram", cex=0.5, font=2)  
#' 
#' @export 
plotTreeDiff <- function(tr1,tr2,tipDiff=NULL,vec1=NULL,vec2=NULL, sizeOfDifferences=FALSE,
                         tipMatch=TRUE, treesFacing=FALSE,
                         baseCol="grey",col1="peachpuff",col2="red2",colourMethod="ramp",palette=lightseasun,
                         ...) {
  
  l <- length(tr1$tip.label)
  
  if (tipMatch) {
    tr2 <- rotateConstr(tr2, tr1$tip.label)
  }
  
  if (is.null(tipDiff)){
    # call tipDiff
    tipDiff <- tipDiff(tr1,tr2,vec1,vec2, sizeOfDifferences)
  }
  
  if (!sizeOfDifferences) {
    # find the number of times each tip appears, in the order that the tip labels are read
    tipSignificance1 <- sapply(tr1$tip.label, function(x)
      tipDiff[which(tipDiff[,1]==x),2])
    tipSignificance2 <- sapply(tr2$tip.label, function(x)
      tipDiff[which(tipDiff[,1]==x),2])
  }
  else {
    # find the size of the tip differences, in the order that the tip labels are read
    tipSignificance1 <- sapply(tr1$tip.label, function(x)
      tipDiff[which(tipDiff[,1]==x),3])
    tipSignificance2 <- sapply(tr2$tip.label, function(x)
      tipDiff[which(tipDiff[,1]==x),3])
  }
  
  # combine:
  significances <- c(tipSignificance1,tipSignificance2)
  
  if (colourMethod=="ramp") {
    colfunc <- colorRampPalette(c(col1,col2))
    if (min(significances)==0) { # make sure tips with no differences are coloured baseCol
      tipSignificance1 <- tipSignificance1 + 1
      tipSignificance2 <- tipSignificance2 + 1
      newSignificances <- c(tipSignificance1,tipSignificance2)
      if (max(significances)==0) {numCols <- 0}
      else {numCols <- max(significances) - min(newSignificances) + 1}
      pal <- c(baseCol,colfunc(numCols))
    }
    else {
      numCols <- max(significances) - min(significances) + 1
      pal <- colfunc(numCols)
    }
  
    tipCols1 <- pal[as.factor(tipSignificance1)]
    tipCols2 <- pal[as.factor(tipSignificance2)]
  }
  
  else {
    tipCols1 <- num2col(tipSignificance1, col.pal=palette)
    tipCols2 <- num2col(tipSignificance2, col.pal=palette)
      if (min(significances)==0) { # make sure tips with no differences are coloured baseCol
        tipCols1[which(tipSignificance1==0)] <- baseCol
        tipCols2[which(tipSignificance2==0)] <- baseCol
      }
  }
    
  if (treesFacing) dir <- "leftwards"
  else dir <- "rightwards"
  
  # plot
  layout(matrix(1:2, 1, 2))
  plot.phylo(tr1, tip.color=tipCols1, no.margin=TRUE, ...)
  plot.phylo(tr2, tip.color=tipCols2, no.margin=TRUE, direction=dir, ...)
  
}

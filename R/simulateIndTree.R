#' Simulate randomised "individuals" tree
#'
#' This function takes in a "category" tree and outputs a simulated corresponding "individuals" tree, for testing the concordance measure.
#'
#' @param catTree object of class phylo, the category-level tree
#' @param itips number of individual tips to assign per category
#' @param permuteCat logical specifying whether to permute the category labels on the category tree before grafting on individual tips. Defaults to FALSE.
#' @param permuteTips logical specifying whether to permute the individual tip labels after building the individual level tree based on the category tree. Defaults to TRUE.
#' @param tipPercent number specifying the percentage of tips to be permuted. Defaults to 100, ignored if permuteTips=FALSE.
#'
#' @import ape
#' @importFrom phytools paste.tree
#' 
#' @examples
#' tree <- simulateIndTree(rtree(3))
#' plot(tree)
#' 
#' @seealso \code{\link{treeConcordance}} \code{\link{makeCollapsedTree}}
#'
#' @export
simulateIndTree <- function(catTree,itips=5,permuteCat=FALSE,permuteTips=TRUE,tipPercent=100) {
  
  if(!is(catTree,"phylo")) stop("The category-level tree catTree must be of class phylo")
  
  bigTree <- catTree # copy of catTree, to have smaller trees grafted on to it
  if (permuteCat){
    bigTree$tip.label <- sample(bigTree$tip.label) # permute tips on category tree, for more discordance
  }
  for (t in catTree$tip.label) {
    bigTree$tip.label[[which(bigTree$tip.label==t)]] <- "NA" # paste.tree will graft a clade into the receptor tree at the "sticky tip" labeled with "NA"
    toglue <- rtree(itips)
    toglue$tip.label <- paste0(t,"_",1:itips)
    toglue$root.edge <- 0
    bigTree <- paste.tree(tr1=bigTree,tr2=toglue)
  }
  if (permuteTips) {
    N <- length(bigTree$tip.label)
    if ((tipPercent<0)||(tipPercent>100)) stop("tipPercent should be a number between 0 and 100")
    p <- floor(N*tipPercent/100)
    samp <- sample(1:N,p)
    tipsToSwap <- bigTree$tip.label[samp]
    bigTree$tip.label[samp] <- sample(tipsToSwap)
  }
  return(bigTree)
}

#' Collapse a tree into a single tip per category
#'
#' Reduce a tree with many tips into a tree with a single tip per category. 
#' Where a category's tips form a monophyletic clade, the clade is replaced by a single tip labelled by that category.
#' Where a category's tips are paraphyletic, the largest clade for that category is treated as above, and all other tips pruned.
#'
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom phangorn Descendants
#'
#' @param tree an object of the class \code{phylo}: the tree to collapse.
#' @param df a two-column data frame linking tip labels (column 2) with their corresponding categories (column 1).
#' @param warnings a logical determining whether a warning should be given if there are paraphyletic categories (default TRUE)
#' 
#' @return
#' A tree (class \code{phylo}) whose tip labels are exactly the set of unique categories from \code{df}.
#' 
#' @seealso \code{\link{treeConcordance}} \code{\link{simulateIndTree}}
#' 
#' @examples  
#' # simulate a tree which is monophyletic per category
#' tree <- simulateIndTree(rtree(5), permuteTips=FALSE)
#' 
#' df <- cbind(sort(rep(rtree(5)$tip.label,5)),sort(tree$tip.label))
#' palette <- c("red","blue","black","green","purple")#' 
#' tipCols <- palette[as.factor(sapply(tree$tip.label, function(l) df[which(df[,2]==l),1]))]
#' plot(tree, tip.color=tipCols)
#' collapsedTree <- makeCollapsedTree(tree,df)
#' plot(collapsedTree, tip.color=palette[as.factor(collapsedTree$tip.label)])
#' 
#' # simulate a tree which is paraphyletic per category
#' tree <- simulateIndTree(rtree(5), tipPercent=20)
#' tipCols <- palette[as.factor(sapply(tree$tip.label, function(l) df[which(df[,2]==l),1]))]
#' plot(tree, tip.color=tipCols)
#' collapsedTree <- makeCollapsedTree(tree,df)
#' plot(collapsedTree, tip.color=palette[as.factor(collapsedTree$tip.label)])
#' 
#' @export 
makeCollapsedTree <- function(tree,df,warnings=TRUE){
  # check whether tips are monophyletic per category,
  # if not, pick a "representative" tip for each category from the largest clade

  categories <- unique(df[,1])
  l <- length(tree$tip.label)
  n <- Nnode(tree)
  tipDescendants <- Descendants(tree, type="tips")
  
  inds_by_category <- split(df[,2], factor(df[,1], categories))
  
  category_is_mono <- sapply(inds_by_category, function(i)
    if (length(i)==1) return(TRUE)
    else any(sapply(tipDescendants, function(j) setequal(i,tree$tip.label[j]) ))
    )

  mono_categories <- categories[category_is_mono]

  tipsToKeep <- vector()
  
  if (any(!category_is_mono)) {
    if(warnings==TRUE)  print("Note: the tree was not monophyletic per category")
    for (c in categories[!category_is_mono]){
      tips <- inds_by_category[[c]]
      # find clade sizes of nodes whose descendant tips all belong to this category 
      clade_sizes_of_this_category <- sapply(1:(n+l), function(x) {
        tmp <- length(intersect(tree$tip.label[tipDescendants[[x]]],tips)) # number of descendants of x which are tips from this category
        if (tmp==length(tipDescendants[[x]])) return(tmp) # if all tip descendants are from this category, return the number of them
        else return(0)
      })
      # keep one tip from the largest clade:
      tipsToKeep <- union(tipsToKeep,tree$tip.label[tipDescendants[[which.max(clade_sizes_of_this_category)[[1]]]]][[1]])
    }
  }
    
    
  tipsToKeep <- union(tipsToKeep,
    sapply(inds_by_category[category_is_mono], function(t) t[[1]]))

  tipsToDrop <- setdiff(tree$tip.label, tipsToKeep)
  collapsedTree <- drop.tip(tree, tip=tipsToDrop, collapse.singles = TRUE) # pruning step

  collapsedTree$tip.label <- as.character(sapply(collapsedTree$tip.label, function(t) df[which(df[,2]==t)[[1]],1])) # relabel the tips by category name

  return(collapsedTree)
}


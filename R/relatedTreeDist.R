#' Tree distance when trees have "related" tips
#'
#' This function calculates the distances between trees whose tips belong to the same categories but are not necessarily identically labelled
#'
#' @param trees an object of class multiphylo
#' @param df a data frame specifying to which category each individual (from all the trees) belongs. Each row gives: an individual (column 2) and its corresponding category (column 1)
#' @param checkTrees a logical (default TRUE) specifying whether the trees should be checked. When TRUE, error messages will be helpful in locating problematic trees, that is, any trees with repeated tip labels, or any trees with missing categories.
#'
#' @import ape
#' @importFrom combinat combn
#' @importFrom fields rdist
#'
#' @examples
#' # we will simulate some trees as an example, each "based" on the same tree:
#' baseTree <- rtree(5)
#' baseTree$tip.label <- letters[5:1]
#' plot(baseTree)
#'
#' tree1 <- simulateIndTree(baseTree, itips=3, permuteTips=FALSE)
#' tree2 <- simulateIndTree(baseTree, itips=4, permuteTips=FALSE)
#' tree3 <- simulateIndTree(baseTree, itips=4, permuteTips=TRUE, tipPercent=20)
#' tree4 <- simulateIndTree(baseTree, itips=4, permuteTips=TRUE, tipPercent=60)
#' tree5 <- simulateIndTree(baseTree, itips=4, permuteTips=TRUE, tipPercent=100)
#' # combine:
#' trees <- list(tree1,tree2,tree3,tree4,tree5)
#'
#' df <- cbind(sort(rep(letters[1:5],4)),sort(paste0(letters[1:5],"_",rep(1:4,5))))
#' head(df)
#'
#' # Find distances:
#' relatedTreeDist(trees,df)
#'
#' # Note that trees 1 and 2 have different numbers of tips but the relationships between those tips
#' # are identical at the category level, hence the related tree distance is 0.
#' # We can see that the distances between trees increase the more the trees are permuted.
#'
#'
#' @export
relatedTreeDist <- function(trees,df,checkTrees=TRUE){
  l <- length(trees)
  # find numbers of leaves and leaf pairs for each tree.
  # If checkTrees=TRUE, check that each tree has no repeated tip labels.
  tree_leaves <- lapply(1:l, function(x) {
      tip_labels <- trees[[x]]$tip.label
      if ((checkTrees=TRUE)&&(anyDuplicated(tip_labels))) stop(paste0("Each tip in tree ",x," must have a unique label"))
      else tip_labels
    })
  tree_leaf_nums <- sapply(tree_leaves, function(x) length(x))

  # to avoid re-computation, make sure each row of df is unique
  df <- unique(df)

  categories <- unique(df[,1]) # category names

  # create a list of vectors: categories with their individuals. Essentially another, easier to parse, representation of df
  inds_by_category <- split(df[,2], factor(df[,1], categories))

  if (checkTrees==TRUE) {
    # check each category is represented at least once in each tree:
    for (c in 1:length(categories)) {
      inds <- inds_by_category[[c]]
      for (i in 1:l) {
        if(length(intersect(inds,tree_leaves[[i]]))==0) stop(paste0("Tree ",i," is missing any individuals from category ",categories[[c]]))
      }
    }
  }

  # find the tip-tip MRCA depths of all individuals, in each tree
  allTreesTipMRCADepths <- lapply(trees, function(x) tipsMRCAdepths(x))

  # get the pairs of categories i,j, just by number
  cat_pair_nums <-  t(combn(1:length(categories),2, fun=c(), simplify=TRUE))

  # for each tree k and each pair of categories i,j,
  # vk_{i,j} = the average depth of the corresponding leaf pairs in tree_k

  vectors <- sapply(1:l, function(t) {  # vector for each tree

    tree_vec <- vector()
    leaf_pairs <- choose(tree_leaf_nums[[t]],2) # number of leaf pairs in this tree

    for (x in 1:length(cat_pair_nums[,1])) { # running through the pairs of categories:
      cati_num <- cat_pair_nums[x,1]
      catj_num <- cat_pair_nums[x,2]

      # list all the pairs of individuals, one from each category, WHICH OCCUR IN THIS TREE
      # (otherwise there's a potential for a huge slow-down)
      ind_pairs <- expand.grid(
          intersect(inds_by_category[[cati_num]],tree_leaves[[t]]),
          intersect(inds_by_category[[catj_num]],tree_leaves[[t]]))

      # we are still fixing x, i.e. the pair of categories,
      # and now we're going to look through all the pairs of individuals from those categories in this tree
      count <- 1
      depths <- vector()
      for (y in 1:length(ind_pairs[,1])) {
        # look for matches to first individual; look through those for matches to second
        # (note, row.match, intersect were slower)
        possible_rows <- c(which(allTreesTipMRCADepths[[t]][,1]==ind_pairs[y,1]),
                           which(allTreesTipMRCADepths[[t]][,2]==ind_pairs[y,1])) # rows which include the first individual
        ijrow <- possible_rows[[c(which(allTreesTipMRCADepths[[t]][possible_rows,1]==ind_pairs[y,2]),
                   which(allTreesTipMRCADepths[[t]][possible_rows,2]==ind_pairs[y,2]))]] # row which includes the first and second individual

        depths[[count]] <- as.integer(allTreesTipMRCADepths[[t]][ijrow,3])
        count <- count + 1
      }
      # for tree t and category pair x, vector[[t]][[x]] is the mean depth of tip pairs in tree t belonging to category pair x
      tree_vec[[x]] <- mean(depths)
    }
    return(tree_vec)
  })

  dist <- rdist(t(data.frame(vectors)))
  return(as.dist(dist))
}

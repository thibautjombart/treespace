#' Tree concordance
#'
#' This function calculates the concordance between a category tree and an individuals tree
#'
#' @param catTree object of class phylo
#' @param indTree object of class phylo
#' @param df data frame specifying to which category each individual belongs. Each row gives an individual (column 2) and its corresponding category (column 1)
#'
#' @import ape
#'
#' @examples
#' # create an example category tree
#' catTree <- read.tree(text="(C:1,(B:1,A:1):1);")
#' plot(catTree)
#'
#' # make individuals tree with complete concordance:
#' indTree1 <- read.tree(text="(((c4,c3),(c2,c1)),((b1,b2),((a3,a2),a1)));")
#' plot(indTree1)
#' 
#' # create data frame linking categories with individuals
#' df <- cbind(c(rep("A",3),rep("B",2),rep("C",4)),sort(indTree1$tip.label))
#' 
#' treeConcordance(catTree,indTree1,df)
#'
#' # make a less concordant tree:
#' indTree2 <- read.tree(text="((((c4,c3),(c2,c1)),b2),(b1,((a3,a2),a1)));")
#' plot(indTree2)
#' treeConcordance(catTree,indTree2,df)
#' 
#' # simulate larger example:
#' catTree <- rtree(10)
#' indTree3 <- simulateIndTree(catTree, tipPercent=10)
#' df <- cbind(sort(rep(catTree$tip.label,5)),sort(indTree3$tip.label))
#' plot(indTree3)
#' treeConcordance(catTree,indTree3,df)
#' 
#' @seealso \code{\link{simulateIndTree}} \code{\link{makeCollapsedTree}}
#'
#' @export
treeConcordance <- function(catTree,indTree,df) {
    
  if(!is(catTree,"phylo")) stop("The category-level tree catTree must be of class phylo")
  if(!is(indTree,"phylo")) stop("The individual-level tree indTree must be of class phylo")
    
  # find the number of pairs of different individuals:
  catTree_leaves <- length(catTree$tip.label)
  indTree_leaves <- length(indTree$tip.label)
  if (length(unique(catTree$tip.label))!=catTree_leaves) stop("Each tip in the category tree must have a unique label")
  if (length(unique(indTree$tip.label))!=indTree_leaves) stop("Each tip in the individuals tree must have a unique label")
  catTree_pairs <- catTree_leaves * (catTree_leaves-1)/2
  indTree_pairs <- indTree_leaves * (indTree_leaves -1)/2

  # avoid unnecessary computation by reducing df to unique rows
  df <- unique(df)
  
  # function to find category corresponding to individual:
  sigma <- function(g) {
    s_row <- which(df[,2]==g)
    s <- as.character(df[s_row,1])
    if (length(s)!= 1) stop("Problem with matching individual to category: check data frame")
    return(s)
  }

  # check category tree has exactly one tip for every category in df,
  # and that each category is represented at least once in indTree:
  categories <- unique(df[,1]) # category names

  # create a list of vectors: categories with their individuals. Essentially another, easier to parse, representation of df
  inds_by_category <- split(df[,2], factor(df[,1], categories))

  if (length(categories)!=catTree_leaves) stop("The number of categories in df does not match the number of leaves in catTree")
  for (c in categories) {
    if (!(c %in% catTree$tip.label)) stop(paste0("Category ",c," from is missing from the catTree"))
    if (length(intersect(inds_by_category[[c]],indTree$tip.label))==0) stop(paste0("indTree is missing any individuals from category ",c))
  }

  # get the tip MRCA depths for both trees
  catTreeDepths <- tipsMRCAdepths(catTree)
  indTreeDepths <- tipsMRCAdepths(indTree)

  # now go through the individuals tree vector, for each column:
  # look up corresponding categories in df using sigma
  # find root distance in category tree vector
  # record whether they match or not

  paircount=0
  con <- vector()
  for (i in 1:indTree_pairs){
    # get the names of the two tips
    indTree_tip1 <- indTreeDepths[i,1]
    indTree_tip2 <- indTreeDepths[i,2]

    # find their category-level labels
    catTree_tips <- c(sigma(indTree_tip1), sigma(indTree_tip2))

    # only count the tips with different category-level labels
    catTree_are_diff <- catTree_tips[1]!=catTree_tips[2]

    if (catTree_are_diff) {
      catTree_dist <- catTreeDepths[
        intersect(which(catTreeDepths[,1] %in% catTree_tips),
                  which(catTreeDepths[,2] %in% catTree_tips))
        ,3]
      indTree_dist <- indTreeDepths[i,3]
      con[[i]] <- indTree_dist==catTree_dist
      paircount <- paircount + 1
    }
    else {
      con[[i]] <- 0
    }
  }
  return(sum(con)/paircount)
}

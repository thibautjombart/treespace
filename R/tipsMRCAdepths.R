#' Tip-tip MRCA depths
#'
#' This function creates a matrix where columns 1 and 2 correspond to tip labels and column 3 gives the depth of the MRCA of that pair of tips.
#' It is strongly based on \code{treeVec} and is used by \code{relatedTreeDist} and \code{treeConcordance} where tip labels belong to "categories".
#'
#' @param tree An object of class phylo
#'
#' @import ape
#' @importFrom combinat combn
#' @importFrom Rcpp evalCpp
#'
#' @examples
#' tree <- rtree(10)
#' plot(tree)
#' tipsMRCAdepths(tree)
#'
#' @export
tipsMRCAdepths <- function(tree) {
  if(class(tree)!="phylo") stop("Tree should be of class phylo")
  if(!is.rooted(tree)) stop("Function is for rooted trees only")
    root <- min(tree$edge[,1])
    if (length(split(tree$edge[,2], factor(tree$edge[,1], root))[[1]]) < 2) stop("Tree root must have at least two descendants")
  tree$edge.length <- rep(1,length(tree$edge))


  num_leaves <- length(tree$tip.label)
  num_edges <- nrow(tree$edge)

  # We work with ordered labels, using this vector to transform indices.
  tip_order <- match(1:num_leaves, order(tree$tip.label))

  # Ordering the edges by first column places the root at the bottom.
  # Descendants will be placed always before parents.
  edge_order <- ape::reorder.phylo(tree, "postorder", index.only=T)
  edges <- tree$edge[edge_order,]
  edge_lengths <- tree$edge.length[edge_order]

  # We annotated the nodes of the tree in this list. In two passes we are going to
  # compute the partition each node induces in the tips (bottom-up pass) and the distance
  # (in branch length and number of branches) from the root to each node (top-down pass).
  annotated_nodes <- list()

  # Bottom up (compute partitions, we store the branch lengths to compute distances
  # to the root on the way down).
  for(i in 1:num_edges) {

    parent <- edges[i,1]
    child <- edges[i,2]

    # Initialization (leaves).
    if(child <= num_leaves) {
      # We translate the index for the sorted labels.
      child <- tip_order[child]
      # Leaves have as children themselves.
      annotated_nodes[[child]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(child))
    }

    # Aggregate the children partitions (only if we are not visiting a leaf).
    aggregated_partitions <- annotated_nodes[[child]]$partitions[[1]]

    if((child > num_leaves)) {
      for(p in 2:length(annotated_nodes[[child]]$partitions))
        aggregated_partitions <- c(aggregated_partitions, annotated_nodes[[child]]$partitions[[p]])
    }

    # Update the branch length on the child.
    annotated_nodes[[child]]$root_distance <- edge_lengths[i]

    # We have not visited this internal node before.
    if(parent > length(annotated_nodes) || is.null(annotated_nodes[[parent]])) {
      # Assume the first time we get the left child partition.
      annotated_nodes[[parent]] <- list(root_distance=NULL, edges_to_root=1, partitions=list(aggregated_partitions))
    }
    # This is not the first time we have visited the node.
    else {
      # We store the next partition of leaves.
      annotated_nodes[[parent]]$partitions[[length(annotated_nodes[[parent]]$partitions)+1]] <- aggregated_partitions
    }
  }

  # Update the distance to the root at the root (i.e. 0)
  # And the number of edges to the root (i.e. 0).
  annotated_nodes[[num_leaves+1]]$root_distance <- 0
  annotated_nodes[[num_leaves+1]]$edges_to_root <- 0

  # Top down, compute distances to the root for each node.
  for(i in num_edges:1) {
    parent <- edges[i,1]
    child <- edges[i,2]

    # If the child is a leaf we translate the index for the sorted labels.
    if(child <= num_leaves)
      child <- tip_order[child]

    annotated_nodes[[child]]$root_distance <- annotated_nodes[[child]]$root_distance + annotated_nodes[[parent]]$root_distance
    annotated_nodes[[child]]$edges_to_root <- annotated_nodes[[child]]$edges_to_root + annotated_nodes[[parent]]$edges_to_root
  }

  # Distance vectors
  vector_length <- (num_leaves*(num_leaves-1)/2) + num_leaves
  length_root_distances <- double(vector_length)
  topological_root_distances <- integer(vector_length)

  # Fill-in the leaves (notice the involved index translation for leaves).
  # Consensus version: final k entries are the heights of leaves, not pendant lengths
  for (i in (vector_length-num_leaves+1):vector_length) {
    topological_root_distances[[i]] <- annotated_nodes[[i-vector_length+num_leaves]]$edges_to_root
  }
  length_root_distances[(vector_length-num_leaves+1):vector_length] <- edge_lengths[match(1:num_leaves, edges[,2])][order(tree$tip.label)]

  # Instead of computing the lexicographic order for the combination pairs assume we
  # are filling in a symmetric distance matrix (using only the triangular upper part).
  # We just need to "roll" the matrix indices into the vector indices.
  # Examples for (k=5)
  # The combination c(1,4) would be located at position 3 on the vector.
  # The combination c(2,1) would be located at position 1 on the vector because d(2,1) = d(1,2).
  # The combination c(2,3) would be located at position 5 on the vector.

  index_offsets <- c(0, cumsum((num_leaves-1):1))

  # This is the slow part, we compute both vectors as gain would be marginal.
  sapply(annotated_nodes, function(node) {

    # We skip leaves and the root (if the MRCA for M groups of leaves is at the root
    # all combinations of leaves -among different groups- have 0 as distance to the root).
    # For large trees this can spare us of computing a lot of combinations.
    # Example: In a perfectly balanced binary tree (N/2 leaves at each side of the root),
    # at the root we'd save (N/2) * (N/2) combinations to update. Worst case scenario is completely
    # unbalanced tree (N-1,1), we'd save in that case only N-1 combinations.

    # Make sure we are not visiting a leaf or the root.
    if(length(node$partitions) > 1 && node$root_distance > 0) {

      # Update all combinations for pairs of leaves from different groups.
      num_groups <- length(node$partitions)
      for(group_a in 1:(num_groups-1)) {
        for(group_b in (group_a+1):num_groups) {
          updateDistancesWithCombinations(length_root_distances, topological_root_distances, node$partitions[[group_a]],
                                          node$partitions[[group_b]], index_offsets, node$root_distance, node$edges_to_root)
        }
      }

    }
  })

  # make tip label vector
  tips_in_order <- tree$tip.label[order(tree$tip.label)]
  label_pairs <- t(combn(tips_in_order,2, fun=c(), simplify=TRUE))
  #self_pairs <- rbind(tips_in_order,tips_in_order)
  tipsMRCAvec <- cbind(label_pairs,topological_root_distances[1:(num_leaves*(num_leaves-1)/2)])
  colnames(tipsMRCAvec) <- c("tip1","tip2","rootdist")

  return(tipsMRCAvec)
}

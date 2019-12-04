#' Find MRCIs
#'
#' Function to find the most recent common infector (MRCI) matrix from "who infected whom" information.
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param wiw a two-column matrix where the first column gives the infectors and the second column gives the infectees; each row corresponds to a transmission event from an infector to an infectee.
#'
#' @return Returns three objects:
#' \itemize{
#' \item \code{sourceCase}: the number of the node which is the source case, i.e. the common infector of all cases (outputs a warning if there is more than one source case).
#' \item \code{mrcis}: a matrix where, for each pair of individuals i and j, the entry (i,j) is the node number of their MRCI. Note that if i infected j then this entry is i itself.
#' \item \code{mrciDepths}: a matrix where, for each pair of individuals i and j, the entry (i,j) is the depth of their MRCI, defined as the number of edges from the source case. The source case has depth zero, its direct infectees have depth 1, and so on. 
#' }
#' 
#' @importFrom compiler cmpfun
#' @importFrom combinat combn2
#'
#' @examples
#'
#' ## a simple who infected whom matrix:
#' tree1 <- cbind(Infector=1:5,Infectee=2:6) 
#' findMRCIs(tree1)
#'
#'
#' @export
findMRCIs <- compiler::cmpfun(function(wiw) { 
  # expect wiw to be a "who infected whom" matrix: column 1 is infectors, column 2 infectees
  if (!is.matrix(wiw)) stop("The who infected whom information supplied should be of class matrix.")
  
  # convert whatever the wiw entries are into integers, but preserve the names for later
  wiwCopy <- matrix(0,nrow(wiw),2)
  wiwNames <- unique(c(wiw[,1],wiw[,2]))
  wiwCopy[,1] <- match(wiw[,1],wiwNames)
  wiwCopy[,2] <- match(wiw[,2],wiwNames)
  wiw <- wiwCopy
  
  initial <- min(wiw) # the number of the first case
  # if the first case is "0" we'll get problems - add one to everything
  if (initial==0) {
    wiw <- wiw + 1
    initial <- min(wiw)
    }
  
  l <- max(wiw) # largest node number (note, not necessarily the "last" case)
  
  # make a reference list where entry i gives all the direct descendant(s) of case i
  DirectDesc <- lapply(initial:l, function(x) {
    wiw[which(wiw[,1]==x),2] # direct infectee(s) of x, if any
  })
  
  numDDs <- sapply(DirectDesc, length)

  AllDesc <- lapply(initial:l, function(x) {
    # for each node x:
    DDs <- DirectDesc[[x]] # find direct descendants
    tmpNumDDs <- length(DDs) # store the number of them
    allDs <- rep(NA, length(initial:l)) # start a vector of all descendants of x, including x itself
    allDs[[1]] <- x
    
    while(tmpNumDDs!=0){
      allDs[which(is.na(allDs))[1:tmpNumDDs]] <- DDs
      DDs <- do.call(c, DirectDesc[DDs]) # now find direct descendants of these, as a vector...
      tmpNumDDs <- length(DDs)
      }
    
    return(allDs[!is.na(allDs)])
  })

  M <- matrix(0, nrow=l, ncol=l) # initialise matrix
  
  # traverse nodes from initial case, down the "tree"
  for (tmp in initial:l){
    
    # tmp is the MRCA of all its descendants with itself
    # (if tmp has no descendants, this syntax is still ok)
    for (i in AllDesc[[tmp]]) {
    M[tmp,i] <- M[i,tmp] <- tmp
    }
    
    # get the direct cases infected by tmp. If there are two or more, then tmp is the MRCA of all pairs of cases descending from different DDs
    tmpDDs <- DirectDesc[[tmp]]
    numTmpDDs <- length(tmpDDs)
      
    if (numTmpDDs>1) { # then tmp is the MRCI of descendant one and descendant two, and of each pair comprised of a descendant from one and a descendant from two
      pairs <- combn2(tmpDDs)
      for (j in 1:length(pairs[,1])) {
        M[AllDesc[[pairs[j,1]]],AllDesc[[pairs[j,2]]]] <- M[AllDesc[[pairs[j,2]]],AllDesc[[pairs[j,1]]]]  <- tmp
        }
    } # end if loop
  }# end for tmp loop
  
  diag(M) <- initial:l # we define the diagonal elements of M to be the nodes themselves
  
  # note: if there are multiple source cases then there will still be zero entries in the matrix, which haven't been overwritten since the initialisation step
  # if there is a single source case, it's the one whose entry in AllDescs is initial:l
  if (min(M)==0) {
    warning("No single source case; returning source of longest transmission chain")
    transChainLengths <- sapply(AllDesc, function(x) length(x))
    sourceOfLongestChain <- which(transChainLengths==max(transChainLengths))
    sourceCase <- sourceOfLongestChain
  }
  else {
    sourceCase <- NA
    count <- 1
    while(is.na(sourceCase)){
      isSource <- setequal(AllDesc[[count]],initial:l)
      if (isSource) {sourceCase <- count}
      count <- count + 1
    }
  }
    
  # finally, we want to know the "depth" of each case (distance from source case)
  # and we want a version of M where the entries are the depths of the mrcis, rather than the IDs of the mrcis
  
  depths <- rep(NA, l)
  count <- 0
  depths[[sourceCase]] <- count
  SCdescs <- unlist(DirectDesc[[sourceCase]], FALSE, FALSE)
  while(length(SCdescs)!=0){
    count <- count + 1
    depths[SCdescs] <- count
    SCdescs <- unlist(DirectDesc[SCdescs], FALSE, FALSE)
  }
  
  D <- matrix(ncol=l, nrow=l, depths[M]) # create matrix D of depths
  M <- matrix(ncol=l, nrow=l, wiwNames[M]) # convert M into the node names

  # re-associate the node names with the columns and rows of the matrices
  colnames(M) <- rownames(M) <- colnames(D) <- rownames(D) <- wiwNames


  return(list(sourceCase=wiwNames[[sourceCase]],mrcis=M,mrciDepths=D))
})

#' Transmission tree distance
#'
#' Function to find the distance between transmission trees by comparing their MRCI depth matrices; to be precise, by finding the Euclidean distance between the tree vectors, restricted to their sampled node entries.
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param matList a list of matrices, each of which is the output of \code{findMRCIs$mrciDepths}
#' @param sampled a vector of node IDs which corresponds to those nodes which are sampled cases. Default is to treat all nodes as sampled cases.
#'
#' @return Returns a distance matrix, where entry (i,j) is the transmission tree distance between matrices i and j in \code{matList} 
#'
#' @importFrom compiler cmpfun
#' @importFrom fields rdist
#'
#' @examples 
#' # create some simple "who infected whom" scenarios:
#' tree1 <- cbind(Infector=1:5,Infectee=2:6) 
#' tree2 <- cbind(Infector=c(1,5,2,2,3),Infectee=2:6)
#' tree3 <- cbind(Infector=c(2,2,3,4,5),Infectee=c(1,3,4,5,6)) 
#' # create list of the MRCI depth matrices:
#' matList <- lapply(list(tree1,tree2,tree3), function(x) findMRCIs(x)$mrciDepths)
#' 
#' # transmission tree distance, assuming all cases are sampled:
#' wiwTreeDist(matList)
#' # transmission tree distance when cases 1, 2 and 4 are sampled:
#' wiwTreeDist(matList, sampled=c(1,2,4))
#'
#' @export
wiwTreeDist <- compiler::cmpfun(function(matList, sampled=NULL) {
  if (is.null(sampled)) {sampled <- 1:length(matList[[1]][1,])}
  
  matVecs <- lapply(matList, function(x) as.vector(x[sampled,sampled]))
  df <- t(data.frame(matVecs))
  
  ## get pairwise Euclidean distances ##
  D <- as.dist(rdist(df))
  
  return(D)
})

#' Median transmission tree
#'
#' Function to find the median of a list of transmission scenarios
#'
#' @author  Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @param matList a list of matrices, each of which is the output of \code{findMRCIs$mrciDepths}
#' @param sampled a vector of node IDs which corresponds to those nodes which are sampled cases
#' @param weights optional vector of weights to correspond to the entries of matList
#'
#' @return Returns three objects:
#' \itemize{
#' \item \code{centre}: the mean of the matList entries, restricted to the sampled cases
#' \item \code{distances}: for each entry of matList, its distance from \code{centre}
#' \item \code{mindist}: the minimum of \code{distances}
#' \item \code{median}: the number of the median entry of matList, i.e. the one(s) which achieve the \code{mindist} from the \code{centre}.
#' }
#'
#' @importFrom compiler cmpfun
#' @importFrom combinat combn2
#'
#' @examples 
#' # create some simple "who infected whom" scenarios:
#' tree1 <- cbind(Infector=1:5,Infectee=2:6) 
#' tree2 <- cbind(Infector=c(1,5,2,2,3),Infectee=2:6)
#' tree3 <- cbind(Infector=c(2,2,3,4,5),Infectee=c(1,3,4,5,6)) 
#' # create list of the MRCI depth matrices:
#' matList <- lapply(list(tree1,tree2,tree3), function(x) findMRCIs(x)$mrciDepths)
#' 
#' # median tree, assuming all cases are sampled:
#' wiwMedTree(matList)
#' # median tree when cases 1, 2 and 4 are sampled:
#' wiwMedTree(matList, sampled=c(1,2,4))
#'
#' @export
wiwMedTree <- compiler::cmpfun(function(matList, sampled=NULL, weights=NULL){
  l <- length(matList)
  
  if (is.null(sampled)) {sampled <- 1:length(matList[[1]][1,])}
  if (is.null(weights)) {weights <- rep(1,l)}
  
  matVecs <- t(sapply(matList, function(x) as.vector(x[sampled,sampled])))
  
  centre <- weights %*% matVecs/l
  
  ## Distances to the centre.
  distances <- apply(matVecs, 1, function(m){sqrt(sum((m-centre)^2))})
  
  ## Get the indices for the median wiw(s).
  min_distance <- min(distances)
  median_wiw <- which(min_distance == distances)
  
  return(list(centre=centre, distances=distances, mindist=min_distance, median=median_wiw))
})


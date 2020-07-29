#'
#' Phylogenetic tree exploration
#'
#' Compares phylogenetic trees using a choice of metrics / measures, and maps their pairwise distances into a small number of dimensions for easy visualisation and identification of clusters.
#'
#' @param x an object of the class multiPhylo
#' @param method the method for summarising the tree as a vector.
#' Choose from:
#' \itemize{
#' \item \code{treeVec} (default) the Kendall Colijn metric vector (for rooted trees)
#' \item \code{BHV} the Billera, Holmes Vogtmann metric using \code{dist.multiPhylo} from package \code{distory} (for rooted trees)
#' \item \code{KF} the Kuhner Felsenstein metric (branch score distance) using \code{KF.dist} from package \code{phangorn} (considers the trees unrooted)
#' \item \code{RF} the Robinson Foulds metric using \code{RF.dist} from package \code{phangorn} (considers the trees unrooted)
#' \item \code{wRF} the weighted Robinson Foulds metric using \code{wRF.dist} from package \code{phangorn} (considers the trees unrooted)
#' \item \code{nNodes} the Steel & Penny tip-tip path difference metric, (topological, ignoring branch lengths), using \code{path.dist} from package \code{phangorn} (considers the trees unrooted)
#' \item \code{patristic} the Steel & Penny tip-tip path difference metric, using branch lengths, calling \code{path.dist} from package \code{phangorn} (considers the trees unrooted)
#' \item \code{CID} the clustering information difference metric, calling \code{ClusteringInformationDistance()} from package \pkg{TreeDist} (considers the trees unrooted)
#' \item \code{PID} the phylogenetic information difference metric, calling \code{PhylogeneticInformationDistance()} from package \pkg{TreeDist} (considers the trees unrooted)
#' \item \code{MS} the matching splits distance, calling \code{MatchingSplitsDistance()} from package \pkg{TreeDist} (considers the trees unrooted)
#' \item \code{MSID} the matching splits information difference metric, calling \code{ClusteringInformationDistance()} from package \pkg{TreeDist} (considers the trees unrooted)
#' \item \code{Abouheif}: performs Abouheif's test, inherited from \code{distTips} in \code{adephylo}. See Pavoine et al. (2008) and \code{adephylo}.
#' \item \code{sumDD}: sum of direct descendants of all nodes on the path, related to Abouheif's test, inherited from \code{distTips} in \code{adephylo}.
#' }
#' @param nf the number of principal components to retain
#' @param lambda a number in [0,1] which specifies the extent to which topology (default, with lambda=0)  or branch lengths (lambda=1) are emphasised in the Kendall Colijn metric.
#' @param return.tree.vectors if using the Kendall Colijn metric, this option will return the tree vectors as part of the output. Note that this can use a lot of memory so defaults to \code{FALSE}.
#' @param processors value (default 1) to be passed to mcmapply specifying the number of cores to use. Must be 1 on Windows (see \code{mcmapply} for more details).
#' @param ... further arguments to be passed to \code{method}.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'
#' @import ape
#' @importFrom ade4 dudi.pco cailliez is.euclid
#' @importFrom adephylo distTips
#' @importFrom distory dist.multiPhylo
#' @importFrom fields rdist
#' @importFrom phangorn KF.dist
#' @importFrom phangorn path.dist
#' @importFrom phangorn RF.dist
#' @importFrom phangorn wRF.dist
#' @importFrom parallel mcmapply
#' @importFrom TreeDist ClusteringInfoDistance MatchingSplitDistance
#'  MatchingSplitInfoDistance PhylogeneticInfoDistance
#'
#' @examples
#'
#' ## generate list of trees
#' x <- rmtree(10, 20)
#' names(x) <- paste("tree", 1:10, sep = "")
#'
#' ## use treespace
#' res <- treespace(x, nf=3)
#' table.paint(as.matrix(res$D))
#' scatter(res$pco)
#'
#' data(woodmiceTrees)
#' woodmiceDists <- treespace(woodmiceTrees,nf=3)
#' plot(woodmiceDists$pco$li[,1],woodmiceDists$pco$li[,2])
#' woodmicedf <- woodmiceDists$pco$li
#' if(require(ggplot2)){
#' woodmiceplot <- ggplot(woodmicedf, aes(x=A1, y=A2)) # create plot
#' woodmiceplot + geom_density2d(colour="gray80") + # contour lines
#' geom_point(size=6, shape=1, colour="gray50") + # grey edges
#' geom_point(size=6, alpha=0.2, colour="navy") + # transparent blue points
#' xlab("") + ylab("") + theme_bw(base_family="") # remove axis labels and grey background
#' }
#'
#' \dontrun{
#' if(require(rgl)){
#' plot3d(woodmicedf[,1], woodmicedf[,2], woodmicedf[,3], type="s", size=1.5,
#' col="navy", alpha=0.5, xlab="", ylab="", zlab="")
#' }
#' }
#'
#'
#' @export
treespace <- function(x, method="treeVec", nf=NULL, lambda=0, return.tree.vectors=FALSE, processors=1, ...){
  
    ## CHECKS ##
    if(!inherits(x, "multiPhylo")) stop("x should be a multiphylo object")
    num_trees <- length(x) # number of trees
    ## fix potential bug with input of two trees
    if(num_trees<3) {
      stop("treespace expects at least three trees. The function treeDist is suitable for comparing two trees.")
    }

    # check for user supplying invalid options (these gave unhelpful error messages before)
    dots <- list(...)
    if(!is.null(dots$return.lambda.function)) stop("return.lambda.function is not compatible with treespace. Consider using multiDist instead.")
    if(!is.null(dots$save.memory)) stop("save.memory is not compatible with treespace. Consider using multiDist instead.")

    # make name labels well defined
    if(is.null(names(x))) names(x) <- 1:num_trees
    else if(length(unique(names(x)))!=num_trees){
      warning("duplicates detected in tree labels - using generic names")
      names(x) <- 1:num_trees
      }
    lab <- names(x)

    # check all trees have same tip labels
    for (i in 1:num_trees) {
      if (!setequal(x[[i]]$tip.label,x[[1]]$tip.label)) {
        stop(paste0("Tree ",lab[[i]]," has different tip labels from the first tree."))
      }
    }

    ## GET DISTANCES BETWEEN TREES, according to method ##
    ## get summary vectors then compute pairwise distances ##
    if (method=="treeVec") {
      df <- t(mcmapply(treeVec, x, lambda=lambda, MoreArgs=dots, mc.cores=processors))
      ## get pairwise Euclidean distances ##
      D <- as.dist(rdist(df))
    }
    else if(method %in% c("Abouheif","sumDD")){
      df <- t(mcmapply(adephylo::distTips, x, method=method, MoreArgs=dots, mc.cores=processors))
      ## get pairwise Euclidean distances ##
      D <- as.dist(rdist(df))
    }
    else if(method=="patristic"){
      D <- path.dist(x, use.weight=TRUE)
    }
    else if(method=="nNodes"){
      D <- path.dist(x, use.weight=FALSE)
    }
    else if(method=="RF"){
      D <- RF.dist(x)
      ## make the distance Euclidean if it isn't already
      if (!ade4::is.euclid(D))  {
        warning("Distance matrix is not Euclidean; making it Euclidean using ade4::cailliez")
        D <- ade4::cailliez(D, print=FALSE)
      }
    }
    else if(method=="wRF"){
      D <- wRF.dist(x)
      ## make the distance Euclidean if it isn't already
      if (!ade4::is.euclid(D))  {
        warning("Distance matrix is not Euclidean; making it Euclidean using ade4::cailliez")
        D <- ade4::cailliez(D, print=FALSE)
      }
    }
    else if(method=="KF"){
      D <- KF.dist(x)
    }
    else if(method=="BHV"){
      D <- dist.multiPhylo(x)
      ## make the distance Euclidean if it isn't already
      if (!ade4::is.euclid(D))  {
        warning("Distance matrix is not Euclidean; making it Euclidean using ade4::cailliez")
        D <- ade4::cailliez(D, print=FALSE)
      }
    }
    else if (method == 'CID') {
      D <- ClusteringInfoDistance(x)
    }
    else if (method == 'PID') {
      D <- PhylogeneticInfoDistance(x)
    }
    else if (method == 'MS') {
      D <- MatchingSplitDistance(x)
    }
    else if (method == 'MSID') {
      D <- MatchingSplitInfoDistance(x)
    }

    ## restore labels
    attr(D,"Labels") <- lab

    ## perform PCoA/MDS ##
    pco <- dudi.pco(D, scannf=is.null(nf), nf=nf)


    ## BUILD RESULT AND RETURN ##
    if (return.tree.vectors==TRUE) {
      out <- list(D=D, pco=pco, vectors=df)
    }
    else {
      out <- list(D=D, pco=pco)
    }
    return(out)
} # end treespace

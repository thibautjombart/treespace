% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transmissionTrees.R
\name{findMRCIs}
\alias{findMRCIs}
\title{Find MRCIs}
\usage{
findMRCIs(wiw)
}
\arguments{
\item{wiw}{a two-column matrix where the first column gives the infectors and the second column gives the infectees; each row corresponds to a transmission event from an infector to an infectee.}
}
\value{
Returns three objects:
\itemize{
\item \code{sourceCase}: the number of the node which is the source case, i.e. the common infector of all cases (outputs a warning if there is more than one source case).
\item \code{mrcis}: a matrix where, for each pair of individuals i and j, the entry (i,j) is the node number of their MRCI. Note that if i infected j then this entry is i itself.
\item \code{mrciDepths}: a matrix where, for each pair of individuals i and j, the entry (i,j) is the depth of their MRCI, defined as the number of edges from the source case. The source case has depth zero, its direct infectees have depth 1, and so on. 
}
}
\description{
Function to find the most recent common infector (MRCI) matrix from "who infected whom" information.
}
\examples{

## a simple who infected whom matrix:
tree1 <- cbind(Infector=1:5,Infectee=2:6) 
findMRCIs(tree1)


}
\author{
Michelle Kendall \email{michelle.louise.kendall@gmail.com}
}

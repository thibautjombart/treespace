[![Build
status](https://ci.appveyor.com/api/projects/status/klr8khh1ieb26rh4/branch/master?svg=true)](https://ci.appveyor.com/project/thibautjombart/treespace/branch/master)
[![CRAN Status
Badge](http://www.r-pkg.org/badges/version/treespace)](https://cran.r-project.org/package=treespace)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/treespace)](https://cran.r-project.org/package=treespace)

*treespace*: exploration of landscapes of phylogenetic trees
============================================================

*treespace* implements new methods for the exploration and analysis of
distributions of phylogenetic trees for a given set of taxa.

Installing *treespace*
----------------------

To install the development version from github:

    library(devtools)
    install_github("thibautjombart/treespace")

The stable version can be installed from CRAN using:

    install.packages("treespace")

Then, to load the package, use:

    library("treespace")

    ## Loading required package: ape

    ## Loading required package: ade4

    ## Creating a generic function for 'toJSON' from package 'jsonlite' in package 'googleVis'

    ## Registered S3 method overwritten by 'spdep':
    ##   method   from
    ##   plot.mst ape

    ## Registered S3 methods overwritten by 'adegraphics':
    ##   method         from
    ##   biplot.dudi    ade4
    ##   kplot.foucart  ade4
    ##   kplot.mcoa     ade4
    ##   kplot.mfa      ade4
    ##   kplot.pta      ade4
    ##   kplot.sepan    ade4
    ##   kplot.statis   ade4
    ##   scatter.coa    ade4
    ##   scatter.dudi   ade4
    ##   scatter.nipals ade4
    ##   scatter.pco    ade4
    ##   score.acm      ade4
    ##   score.mix      ade4
    ##   score.pca      ade4
    ##   screeplot.dudi ade4

Content overview
----------------

The main functions implemented in *treespace* are:

-   **`treespace`**: explore landscapes of phylogenetic trees

-   **`treespaceServer`**: open up an application in a web browser for
    an interactive exploration of the diversity in a set of trees

-   **`findGroves`**: identify clusters of similar trees

-   **`plotGroves`**: scatterplot of groups of trees, and
    **`plotGrovesD3`** which enables interactive plotting based on d3.js

-   **`medTree`**: find geometric median tree(s) to summarise a group of
    trees

-   **`wiwTreeDist`**: find the distance between transmission trees by
    comparing their MRCI depth matrices

-   **`wiwMedTree`**: find the median of a list of transmission
    scenarios

-   **`relatedTreeDist`**: calculate the distances between trees whose
    tips belong to the same categories but are not necessarily
    identically labelled

-   **`treeConcordance`**: calculate the concordance between a category
    tree and an individuals tree

Other functions are central to the computations of distances between
trees:

-   **`treeVec`**: characterise a tree by a vector

-   **`treeDist`**: find the distance between two tree vectors

-   **`multiDist`**: find the pairwise distances of a list of trees

-   **`refTreeDist`**: find the distances of a list of trees from a
    reference tree

-   **`tipDiff`**: for a pair of trees, list the tips with differing
    ancestry

-   **`plotTreeDiff`**: plot a pair of trees, highlighting the tips with
    differing ancestry

-   **`findMRCIs`**: find the most recent common infector (MRCI) matrix
    from “who infected whom” information

-   **`tipsMRCAdepths`**: similar to `treeVec` but the output is a
    matrix where columns 1 and 2 correspond to tip labels and column 3
    gives the depth of the MRCA of that pair of tips

Distributed datasets include:

-   **`woodmiceTrees`**: illustrative set of 201 trees built using the
    neighbour-joining and bootstrapping example from the *woodmice*
    dataset in the *ape* documentation.

-   **`DengueTrees`**: 500 trees sampled from a BEAST posterior set of
    trees from (Drummond and Rambaut, 2007)

-   **`DengueSeqs`**: 17 dengue virus serotype 4 sequences from
    (Lanciotti *et al.*, 1997), from which the `DengueTrees` were
    inferred.

-   **`DengueBEASTMCC`**: the maximum clade credibility (MCC) tree from
    the `DengueTrees`.

Documentation
-------------

*treespace* comes with the following vignettes:

-   [*introduction*](https://cran.r-project.org/package=treespace/vignettes/introduction.html):
    general introduction using a worked example.

-   [*Dengue
    example*](https://cran.r-project.org/package=treespace/vignettes/DengueVignette.html):
    worked example using a Dengue dataset, used in the *treespace*
    publication.

-   [*transmission
    trees*](https://cran.r-project.org/package=treespace/vignettes/TransmissionTreesVignette.html):
    worked example using transmission trees.

-   [*tip
    categories*](https://cran.r-project.org/package=treespace/vignettes/tipCategories.html):
    introduction to the measures for comparing trees with shared tip
    label “categories”

Contributing / asking a question
--------------------------------

Contributions are welcome via **pull requests**.

Please note that this project is released with a [Contributor Code of
Conduct](https://thibautjombart.github.io/treespace/CONDUCT.html). By
participating in this project you agree to abide by its terms.

Questions, feature requests and bugs can be reported using the package’s
[issue system](https://github.com/thibautjombart/treespace/issues).

Authors
-------

Authors:

-   [Thibaut Jombart](https://thibautjombart.netlify.com/)

-   [Michelle Kendall](https://michellekendall.github.io/)

Contributors:

-   [Jacob
    Almagro-Garcia](http://www.well.ox.ac.uk/jacob-almagro-garcia)

-   [Caroline Colijn](http://www.imperial.ac.uk/people/c.colijn)

Maintainer of the CRAN version:

-   [Michelle Kendall](https://michellekendall.github.io/)

See details of contributions on: <br>
<a href="https://github.com/thibautjombart/treespace/graphs/contributors" class="uri">https://github.com/thibautjombart/treespace/graphs/contributors</a>

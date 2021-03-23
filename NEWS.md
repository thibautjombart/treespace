treespace v1.1.4.1
==================

Fixing a bug which was tripping up our automatic tests: if d is a distance matrix obtained from stats::as.dist() then its i'th element can no longer be accessed with d[[i]] and instead d[i] should be used.

treespace v1.1.4.0
==================

We are very grateful to Joseph Tsui and Alexis Robert from the London School of Hygiene and Tropical Medicine for spotting an error in the way tip labels were handled in calculating transmission tree distances - this has now been fixed and a test added.

treespace v1.1.3.2
==================

This patch corrects the bad usage of "class" which was causing errors in some r-devel-linux versions, and now has a fully specified URL linking from the readme to the code of conduct.

treespace v1.1.3.1
==================

We are very grateful to Kurt Hornik for supplying this fix while Michelle was on maternity leave; the fix handles the update to the R base function `sample`.

treespace v1.1.3
==================

We have updated the dependency details for the package `rgl` to ensure that the treespace Shiny app works on all operating systems: removed a few lines from a vignette which caused warnings on some operating systems; added a reference to the DESCRIPTION file as requested by Uwe, plus and URL and BugReports fields.

treespace v1.1.2
==================

We have updated the maintainer email address and fixed the issue of stating a dependency on an R version with patchlevel non-zero.

treespace v1.1.1
==================

We made minor changes to package dependencies and removed support for parallel computation to satisfy CRAN policies.

treespace v1.1.0
==================

We have added functions for comparing trees with "related" tip sets.

treespace v1.0.1
=================

We made minor improvements to functions and updated some package dependencies. Michelle Kendall became the maintainer of the package again, with thanks to Thibaut Jombart for maintaining it during her leave. 

treespace v1.0.0
================

New package, replacing the previous `treescape`.
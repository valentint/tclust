##  VT::11.10.2023 - this will render the output independent
##  from the version of the package
suppressPackageStartupMessages(library(tclust))

require(tclust)
require(MASS)
#--- EXAMPLE 1 ------------------------------------------

set.seed(123)
sig <- diag (2)
cen <- rep (1,2)
x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
           MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
           MASS::mvrnorm(100, cen * 2.5, sig * 50)
           )

# Two groups and 10% trimming level
(clus <- tclust(x, k = 2, alpha = 0.1, restr.fact = 8))


# Three groups (one of them very scattered) and 0% trimming level
(clus <- tclust(x, k = 3, alpha=0.0, restr.fact = 100))


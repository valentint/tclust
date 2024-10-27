##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'  LG5data data
#'
#' A data set in dimension 10 with three clusters around affine subspaces
#'  of common intrinsic dimension. A 10\% background noise is added uniformly 
#'  distributed in a rectangle containing the three main clusters. 
#' @name LG5data
#' @docType data
#' @usage data(LG5data)
#' @format The first 10 columns are the variables. The last column is the true 
#'  classification vector where symbol "0" stands for the contaminating data points. 
#' @examples 
#' #--- EXAMPLE 1 ------------------------------------------ 
#' data (LG5data)
#' x <- LG5data[, 1:10]
#' clus <- rlg(x, d = c(2,2,2), alpha=0.1, trace=TRUE)
#' plot(x, col=clus$cluster+1)
#' @keywords datasets
NULL

#'  M5data data
#'
#' A bivariate data set obtained from three normal bivariate distributions with 
#'  different scales and proportions 1:2:2. One of the components is very overlapped 
#'  with another one. A 10\% background noise is added uniformly distributed in a rectangle 
#'  containing the three normal components and not very overlapped with the three mixture 
#'  components. A precise description of the M5 data set can be found in 
#'  García-Escudero et al. (2008).
#'
#' @name M5data
#' @docType data
#' @usage data(M5data)
#' @format The first two columns are the two variables. The last column is the true 
#'  classification vector where symbol "0" stands for the contaminating data points. 
#' @source García-Escudero, L.A.; Gordaliza, A.; Matrán, C. and Mayo-Iscar, A. (2008), 
#'  "A General Trimming Approach to Robust Cluster Analysis". Annals of Statistics, 
#'  Vol.36, pp. 1324-1345.
#' @examples 
#' #--- EXAMPLE 1 ------------------------------------------ 
#' data (M5data) 
#' x <- M5data[, 1:2] 
#' clus <- tclust(x, k=3, alpha=0.1, nstart=200, niter1=3, niter2=17, 
#'    nkeep=10, opt="HARD", equal.weights=FALSE, restr.fact=50, trace=TRUE) 
#' plot (x, col=clus$cluster+1)
#' @keywords datasets
NULL

#' Wholesale customers dataset
#'
#' The data set refers to clients of a wholesale distributor. It includes the annual 
#'  spending in monetary units  on diverse product categories.
#'
#' @name wholesale
#' @docType data
#' @usage data(wholesale)
#' @format A data frame containing 440 observations in 8 variables (6 numerical and two categorical). 
#' The variables are as follows:
#'
#' \itemize{
#'      \item \code{Region} Customers' Region - Lisbon (coded as 1), Porto (coded as 2) or Other (coded as 3)
#'      \item \code{Fresh} Annual spending on fresh products
#'      \item \code{Milk} Annual spending on milk products
#'      \item \code{Grocery} Annual spending on grocery products
#'      \item \code{Frozen} Annual spending on frozen products
#'      \item \code{Detergents} Annual spending on detergents and paper products
#'      \item \code{Delicatessen} Annual spending on and delicatessen products
#'      \item \code{Channel}  Customers' Channel - Horeca (Hotel/Restaurant/Café) or 
#'              Retail channel. Horeca is coded as 1 and Retail channel is coded as 2
#'  }
#'
#' @source Abreu, N. (2011). Analise do perfil do cliente Recheio e desenvolvimento de 
#'  um sistema promocional. Mestrado em Marketing, ISCTE-IUL, Lisbon. 
#'  url={https://api.semanticscholar.org/CorpusID:124027622}
#'
#' @examples 
#' #--- EXAMPLE 1 ------------------------------------------ 
#' data (wholesale) 
#' x <- wholesale[, -c(1, ncol(wholesale))] 
#' clus <- tclust(x, k=3, alpha=0.1, nstart=200, niter1=3, niter2=17, 
#'    nkeep=10, opt="HARD", equal.weights=FALSE, restr.fact=50, trace=TRUE) 
#'  plot (x, col=clus$cluster+1)
#'  plot(clus)
#' @keywords datasets
NULL

#' Pinus nigra dataset
#'
#' To study the growth of the wood mass in a cultivated forest of \emph{Pinus nigra} 
#'  located in the north of Palencia (Spain), a sample of 362 trees was studied. 
#'  The data set is made of measurements of heights (in meters), in variable "HT", 
#'  and diameters (in millimetres), in variable "Diameter", of these trees. 
#'  The presence of three linear groups can be guessed apart from a small group 
#'  of trees forming its own cluster with larger heights and diameters one isolated 
#'  tree with the largest diameter but small height. More details on the 
#'  interpretation of this dataset in García-Escudero et al (2010).
#'
#' @name pine
#' @docType data
#' @usage data(pine)
#' @format A data frame containing 362 observations in 2 variables. 
#' The variables are as follows:
#'
#' \itemize{
#'      \item \code{Diameter} Diameter
#'      \item \code{HT} Height
#'  }
#'
#' @references García-Escudero, L. A., Gordaliza, A., Mayo-Iscar, A., and San Martín, R. (2010). 
#'  Robust clusterwise linear regression through trimming. 
#'  \emph{Computational Statistics & Data Analysis}, 54(12), 3057--3069.
#'
#' @keywords datasets
NULL

#' Old Faithful Geyser Data
#'
#' A bivariate data set obtained from the Old Faithful Geyser, containing the 
#'  eruption length and the length of the previous eruption for 271 eruptions 
#'  of this geyser in minutes.
#'
#' @name geyser2
#' @docType data
#' @usage data(geyser2)
#' @format A data frame containing 272 observations in 2 variables. 
#' The variables are as follows:
#'
#' \itemize{
#'      \item \code{Eruption length} The eruption length in minutes. 
#'      \item \code{Previous eruption length} The length of the previous eruption in minutes. 
#'  }
#'
#' @source
#'  This particular data structure can be obtained by applying the following code
#'  to the "Old Faithful Geyser" (\code{faithful} data set (Härdle 1991) in the 
#'  package \code{datasets}):
#'  \cr
#'  \code{f1 <- faithful[,1]}\cr
#'  \code{geyser2 <- cbind(f1[-length(f1)], f1[-1])}\cr
#'  \code{colnames(geyser2) <- c("Eruption length",}\cr
#'  \code{                      "Previous eruption length")}
#
#' @references 
#'  García-Escudero, L.A. and Gordaliza, A. (1999). 
#'  Robustness properties of k-means and trimmed k-means. \emph{Journal of the American Statistical Assoc.}, Vol.94, No.447, 956--969.
#'
#'  Härdle, W. (1991). \emph{Smoothing Techniques with Implementation in S.}, New York: Springer.
#'
#' @keywords datasets
NULL

#' Swiss banknotes data
#'
#' Six variables measured on 100 genuine and 100 counterfeit old Swiss 1000-franc 
#'  bank notes (Flury and Riedwyl, 1988).
#'
#' @name swissbank
#' @docType data
#' @usage data(swissbank)
#' @format A data frame containing 200 observations in 6 variables. 
#' The variables are as follows:
#'
#' \itemize{
#'      \item \code{Length} Length of the bank note 
#'      \item \code{Ht_Left} Height of the bank note, measured on the left 
#'      \item \code{Ht_Right} Height of the bank note, measured on the right 
#'      \item \code{IF_Lower} Distance of inner frame to the lower border 
#'      \item \code{IF_Upper} Distance of inner frame to the upper border 
#'      \item \code{Diagonal} Length of the diagonal 
#'  }
#'
#' @details 
#'  Observations 1--100 are the genuine bank notes and the other 100 observations are the counterfeit bank notes. 
#'
#' @source
#'  Flury, B. and Riedwyl, H. (1988). \emph{Multivariate Statistics, A Practical Approach}, Cambridge University Press. 
#
#' @keywords datasets
NULL
#' Flea
#'
#' Flea-beetle measurements
#'
#' @name flea
#' @docType data
#' @usage data(flea)
#' @format A data frame with 74 rows and 7 variables: six explanatory and one response variable - \code{species}.
#' The variables are as follows:
#'
#' \itemize{
#'   \item tars1: width of the first joint of the first tarsus in microns (the sum of measurements for both tarsi)
#'   \item tars2: the same for the second joint
#'   \item head: the maximal width of the head between the external edges of the eyes in 0.01 mm
#'   \item ade1: the maximal width of the aedeagus in the fore-part in microns
#'   \item ade2: the front angle of the aedeagus ( 1 unit = 7.5 degrees)
#'   \item ade3: the aedeagus width from the side in microns
#'   \item species, which species is being examined - \code{Concinna}, \code{Heptapotamica}, \code{Heikertingeri}
#' }
#'
#'
#' @references A. A. Lubischew (1962), On the Use of Discriminant Functions in Taxonomy, \emph{Biometrics}, \bold{18}4 pp.455--477.
#'
#' @examples
#'  data(flea)
#'  head(flea)
#'
#' @keywords datasets
NULL

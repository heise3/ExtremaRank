\name{golub}
\alias{golub}
\alias{golub.cl}
\alias{golub.gnames}

\title{Gene expression dataset from Golub et al. (1999)}

\usage{
data(golub)
}

\description{
  Gene expression data (3051 genes and 38 tumor mRNA samples) from the
  leukemia microarray study of Golub et al. (1999). Pre-processing
  was done as described in Dudoit et al. (2002). The R code for pre-processing is available
in the file \url{../doc/golub.R}.}

\value{
  \item{golub}{matrix of gene expression levels for the 38 tumor mRNA samples, rows correspond to genes (3051 genes) and columns to mRNA samples.}
  \item{golub.cl}{numeric vector indicating the tumor class, 27 acute lymphoblastic leukemia (ALL) cases (code 0) and 11 acute myeloid leukemia (AML) cases (code 1). }
  \item{golub.gnames}{a matrix containing the names of the 3051 genes for the expression matrix \code{golub}. The three columns correspond to the gene \code{index}, \code{ID}, and \code{Name}, respectively.
  }
}

\source{Golub et al. (1999). Molecular classification of cancer: class
discovery and class prediction by gene expression
monitoring, \emph{Science}, Vol. 286:531-537.\cr
\url{http://www-genome.wi.mit.edu/MPR/}
.}

\references{S. Dudoit, J. Fridlyand, and T. P. Speed (2002). Comparison of discrimination methods for the  classification of tumors using gene expression data. \emph{Journal of the American Statistical Association}, Vol. 97, No. 457, p. 77--87. }
\keyword{datasets}  

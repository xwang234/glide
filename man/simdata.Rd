\name{simdata}
\alias{simdata}
\alias{coeff}
\alias{simdat}
\docType{data}
\title{An example dataset to demostrate the usage of glide
}
\description{
  An example simulation dataset used to demostrate the usage of glide. It includes a dataframe "simdat" and a dataframe "coeff"
}
\usage{data("simdata")}
\format{
coeff is a one-column dataframe. It consists of exposure coefficients of genotypes of 75 snps.
simdat is a data frame consisting of 20000 observations. It has 81 columns, with 75 columns for genotypes of 75 snps and 6 columns of following covariates:
  \describe{
    \item{\code{age}}{an integer vector to indicate age}
    \item{\code{sex}}{a factor vector to indicate gender}
    \item{\code{pc1}}{a numeric vector of the first principal component of genotypes}
    \item{\code{pc2}}{a numeric vector of the second principal component of genotypes}
    \item{\code{pc3}}{a numeric vector of the third principal component of genotypes}
    \item{\code{outcome}}{an binary vector to indicate the outcome, 1/0: the observation has/hasn't cancer}
  }
}
\examples{
  data("simdata")
  ls()
}
\keyword{datasets}
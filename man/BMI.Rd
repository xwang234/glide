\name{BMI}
\alias{BMI}
\docType{data}
\title{An example dataset to demostrate the usage of glide
}
\description{
  An example dataset borrowed from Genetics and Epidemiology of Colorectal Cancer Consortium(GECCO) used to demostrate the usage of glide. It includes a dataframe "dat" and a vector "exposure_coeff"
}
\usage{data("BMI")}
\format{
exposure_coeff is a named vector. It consists of exposure coefficients of genotypes of 77 snps.
dat is a data frame consisting of 20515 observations. It has 84 columns, with 77 columns for genotypes of 77 snps and 7 columns of following covariates:
  \describe{
    \item{\code{gecco_study}}{a factor vector to indicate which study the observation came from}
    \item{\code{age_ref}}{an integer vector to indicate age}
    \item{\code{sex}}{a factor vector to indicate gender}
    \item{\code{pc1}}{a numeric vector of the first principal component of genotypes}
    \item{\code{pc2}}{a numeric vector of the second principal component of genotypes}
    \item{\code{pc3}}{a numeric vector of the third principal component of genotypes}
    \item{\code{outcome}}{an binary vector to indicate the outcome, 1/0: the observation has/hasn't cancer}
  }
}
\examples{
  data("BMI")
  head(exposure_coeff)
  dim(dat)
}
\keyword{datasets}
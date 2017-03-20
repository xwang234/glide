\name{checkdata}
\alias{checkdata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function used to check the core inputs of function glide
}
\description{
It is used to check the core inputs of function glide
}
\usage{
checkdata(formula,exposure_coeff,genotype_columns,data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A an object of class "formula": outcome ~ adusting covariates.
}
  \item{exposure_coeff}{
A named mumeric vector for the external coefficients for genotype-exposure associations.
}
  \item{genotype_columns}{
An integer vector, column indices of data in which genotype data are stored.  
}
  \item{data}{
A data frame A data frame with following components.
}
}


\references{
%% ~put references to the literature/web site here ~
}
\author{
James Y. Dai and X. Wang
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##load an example data; it includes a dataframe "dat" and a named numeric vector  exposure_coeff
data(BMI)
formula=as.formula("outcome~gecco_study+age_ref+sex+pc1+pc2+pc3")
valid=checkdata(formula=formula,exposure_coeff=exposure_coeff,data=dat)
}
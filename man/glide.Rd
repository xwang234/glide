\name{glide}
\alias{glide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to perform global and individual tests for assessing direct effects
}
\description{
A function to perform global and individual tests for assessing direct effects
}
\usage{
glide(formula,exposure_coeff=NULL,genotype_columns=NULL,data,
np=100000,qcutoff=0.2,parallel=TRUE,corenumber=4,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A an object of class "formula": outcome ~ adusting covariates.
}
  \item{exposure_coeff}{
A named mumeric vector for the external coefficients for genotype-exposure associations. 
All the names must be found in "genotype_columns" columns of data.
}
  \item{genotype_columns}{
An integer vector, column indices of data in which genotype data are stored.  
}
  \item{data}{
A data frame A data frame with following components.
}
  \item{np}{
A integer for the number of parametric simulations being generated compute q-values.
}
  \item{qcutoff}{
A integer for the q-value cutoff to declare individual variant’s significance.
}
  \item{parallel}{
If TRUE, use parallel foreach to speed up computation.
}
  \item{corenumber}{
A integer to set the number of computer cores used in parallel computation.
}
  \item{verbose}{
If TRUE, print time-spending info.
}
}
\details{

}
\value{
%%  ~Describe the value returned
A data frame containing observed p-value, expected null p-value,fwer,q-value of each snp.
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
##load an example data; it includes a dataframe "dat" and a named numeric vector  "exposure_coeff"
data(BMI)
formula=as.formula("outcome~as.factor(gecco_study)+age_ref+sex+pc1+pc2+pc3")
genotype_columns=which(grepl("^rs",colnames(dat)))
out=glide(formula=formula,exposure_coeff=exposure_coeff,genotype_columns,data=dat,np=100000,qcutoff=0.2,parallel=TRUE,corenumber=4,verbose=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Causal Analysis,Direct Effects}

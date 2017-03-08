\name{glide}
\alias{glide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to perform global and individual tests for direct effects
}
\description{
A function to perform global and individual tests for direct effects
}
\usage{
glide(data, genotype_columns = NULL, adjusting_covariate_columns = NULL, outcome_column = NULL, coeff = NULL, np = 1e+05, qcutoff = 0.2, parallel = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame used to access the following data.
}
  \item{genotype_columns}{
An integer vector, column indices of data in which genotyp data are stored.  
}
  \item{adjusting_covariate_columns}{
An integer vector, column indices of data in which adjusting covariates are stored.  
}
  \item{outcome_column}{
An integer, in which column of data the outcome is stored. 
}
  \item{coeff}{
A named mumeric vector for the external coefficients.
}
  \item{np}{
A integer for the number of permutation used to compute q-values.
}
  \item{qcutoff}{
A integer for the q-value cutoff.
}
  \item{parallel}{
If TRUE, use parallel foreach to speed up computation. Must register parallel before hand, such as doMC or others. See the example below.
}
  \item{verbose}{
If TRUE, print time-spending info.
}
}
\details{

}
\value{
%%  ~Describe the value returned
A data frame containing observed p-value, expected null p-value,fwer,q-value of each genotype.
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
##load an example data
data(BMI)
dat$gecco_study=as.factor(dat$gecco_study)

genotype_columns=which(colnames(dat) \%in\% names(coeff))
adjusting_covariate_columns=which(colnames(dat) \%in\% c("pc1","pc2","pc3","age_ref","sex","gecco_study"))
outcome_column=which(colnames(dat)=="outcome")
result=glide(data=dat,genotype_columns=genotype_columns,adjusting_covariate_columns=adjusting_covariate_columns,outcome_column=outcome_column,coeff=coeff)
\dontrun{
#Parallel
library(doMC)
registerDoMC(cores=4)
result=glide(data=dat,genotype_columns=genotype_columns,adjusting_covariate_columns=adjusting_covariate_columns,outcome_column=outcome_column,coeff=coeff,parallel=T)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Direct Effects}

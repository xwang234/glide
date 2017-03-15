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
glide(data,genotype_columns=NULL,adjusting_covariate_columns=NULL,
                  outcome_column=NULL,exposure_coeff=NULL,np=100000,qcutoff=0.2,parallel=TRUE,corenumber=4,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
A data frame A data frame with following components.
}
  \item{genotype_columns}{
An integer vector, column indices of data in which genotyp data are stored.  
}
  \item{adjusting_covariate_columns}{
An integer vector, column indices of data in which adjusting covariates are stored.  
}
  \item{outcome_column}{
An integer, in which column of data the case-control outcome is stored. 
}
  \item{exposure_coeff}{
A named mumeric vector for the external coefficients for genotype-exposure associations.
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
A integer to set the number of computer cores participating in parallel computation.
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
genotype_columns=which(colnames(dat) \%in\% names(exposure_coeff))
adjusting_covariate_columns=which(colnames(dat) \%in\% c("pc1","pc2","pc3","age_ref","sex","gecco_study"))
outcome_column=which(colnames(dat)=="outcome")
out=glide(data=dat,genotype_columns=genotype_columns,adjusting_covariate_columns=adjusting_covariate_columns,outcome_column=outcome_column,exposure_coeff=exposure_coeff)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Causal Analysis,Direct Effects}

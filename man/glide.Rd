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
np=100000,qcutoff=0.2,parallel=TRUE,corenumber=1,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A an object of class "formula": outcome ~ adusting covariates.
}
  \item{exposure_coeff}{
A named mumeric vector or a one-column dataframe for the external coefficients 
for genotype-exposure associations. 
All the names/rownames must be found in "genotype_columns" columns of data.
}
  \item{genotype_columns}{
An integer vector, column indices of data in which genotype data are stored.  
}
  \item{data}{
A data frame composed of formula-related components and snp genotype columns.
}
  \item{np}{
A integer for the number of parametric simulations being generated compute q-values.
}
  \item{qcutoff}{
A integer for the q-value cutoff to declare significance of individual variants.
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
A data frame containing observed p-value, expected null p-value,fwer,q-value of each snp.
}
\value{
%%  ~Describe the value returned
A data frame containing observed p-value, expected null p-value,fwer,q-value of each snp.
}
\references{
J. Y. Dai, U. Peters, X. Wang, J. Kocarnik, J. et al. Diagnostics of pleiotropy in mendelian
randomization studies: Global and individual tests for direct effects. The American
Journal of Human Genetics, in press.

}
\author{
James Y. Dai and X. Wang
}

\examples{
##load an example data; it includes a dataframe "simdat" and a one-column dataframe  "coeff"
data(simdata)
formula=as.formula("outcome~age+sex+pc1+pc2+pc3")
genotype_columns=which(grepl("^SNP",colnames(simdat)))
out=glide(formula=formula,exposure_coeff=coeff,genotype_columns,data=simdat,np=100000,
qcutoff=0.2,parallel=TRUE,corenumber=1,verbose=TRUE)
##use more cores
\dontrun{
out=glide(formula=formula,exposure_coeff=coeff,genotype_columns,data=simdat,np=100000,
qcutoff=0.2,parallel=TRUE,corenumber=4,verbose=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Causal Inference, Mendelian Randomization, Pleiotropy}

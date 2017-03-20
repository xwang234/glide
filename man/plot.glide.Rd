\name{plot.glide}
\alias{plot.glide}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function used to draw the q-q plot 
}
\description{
It is used to draw the q-q plot
}
\usage{
glide_plot(out,qcutoff=0.2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{out}{
out is a dataframe, it is the outcome of the glide function.
}
  \item{qcutoff}{
qcutoff is the q-value cutoff.
}
}

\author{
James Y. Dai and X. Wang
}

\examples{
#first generate output using the glide function
##load an example data; it includes a dataframe "dat" and a named numeric vector  exposure_coeff
data(BMI)
formula=as.formula("outcome~as.factor(gecco_study)+age_ref+sex+pc1+pc2+pc3")
genotype_columns=which(grepl("^rs",colnames(dat)))
out=glide(formula=formula,exposure_coeff=exposure_coeff,genotype_columns,data=dat,np=100000,qcutoff=0.2,parallel=TRUE,corenumber=4,verbose=TRUE)

##draw the plot
plot.glide(out)
}
\name{plot.egger}
\alias{plot.egger}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function used to draw Egger plot
}
\description{
It is used to draw MR-Egger regression plot.
}
\usage{
\method{plot}{egger}(x,qcutoff=0.2,xlab="Genetic association with the exposure",
                       ylab="Genetic association with the outcome",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
x is a dataframe, it is the outcome of the glide function.
}
\item{qcutoff}{
qcutoff is the q-value cutoff.
}
\item{xlab}{
xlab is the default x label.
}
\item{ylab}{
ylab is the default y label.
}
\item{...}{
the rest of plot arguments.
}
}

\author{
James Y. Dai and X. Wang
}

\examples{
\dontrun{
##first generate output using the glide function
#load an example data
data(simdata)
simdat=simdata$simdat
coeff=simdata$coeff

formula=as.formula("outcome~age+sex+pc1+pc2+pc3")
genotype_columns=which(grepl("^SNP",colnames(simdat)))

#run glide
out=glide(formula=formula,exposure_coeff=coeff,genotype_columns,data=simdat,np=100000,
qcutoff=0.2,parallel=TRUE,corenumber=1,verbose=TRUE)

##draw the Egger plot
plot.egger(out)
}
}
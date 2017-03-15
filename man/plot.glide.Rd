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
##load an example data
data(BMI)
dat$gecco_study=as.factor(dat$gecco_study)
genotype_columns=which(colnames(dat) \%in\% names(exposure_coeff))
adjusting_covariate_columns=which(colnames(dat) \%in\% c("pc1","pc2","pc3","age_ref","sex","gecco_study"))
outcome_column=which(colnames(dat)=="outcome")
out=glide(data=dat,genotype_columns=genotype_columns,adjusting_covariate_columns=adjusting_covariate_columns,outcome_column=outcome_column,exposure_coeff=exposure_coeff)
##draw the plot
plot.glide(out)
}
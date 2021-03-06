% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIT.R
\name{init}
\alias{init}
\title{A raw API for initializing model parameters.}
\usage{
init(expression, weight, attribute, weather, recipe)
}
\arguments{
\item{expression}{An object that represents gene expression data.
The object can be created from a dumped/saved dataframe
of size \code{nsamples * ngenes}
using \code{FIT::load.expression()}.
(At the moment it is an instance of a hidden class IO$Attribute,
but this may be subject to change.)}

\item{weight}{A matrix of size \code{nsamples * ngenes}
    that during regression penalizes errors from each sample 
    using the formula
    \code{sum_{s in samples} (weight_s) (error_s)^2}.

    Note that, unlike for \code{FIT::train()}, this argument
    is NOT optional.}

\item{attribute}{An object that represents the attributes of a
microarray/RNA-seq data.
The object can be created from a dumped/saved dataframe
of size \code{nsamples * nattributes}
using \code{FIT::load.attribute()}.
(At the moment it is an instance of a hidden class IO$Attribute,
but this may be subject to change.)}

\item{weather}{An object that represents actual or hypothetical weather data
with which the training of models are done.
The object can be created from a dumped/saved dataframe
of size \code{ntimepoints * nfactors}
using \code{FIT::load.weather()}.
(At the moment it is an instance of a hidden class IO$Weather,
but this may be subject to change.)}

\item{recipe}{An object that represents the training protocol of models.
A recipe can be created using \code{FIT::make.recipe()}.}
}
\value{
A collection of models whose parameters are
    set by using the \code{'init'} method in the argument \code{recipe}.
}
\description{
Note: use \code{train()} unless the user is willing to
accept breaking API changes in the future.
}

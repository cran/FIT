% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FIT.R
\name{convert.attribute}
\alias{convert.attribute}
\title{Converts attribute data from a dataframe into an object.}
\usage{
convert.attribute(data, sample = NULL)
}
\arguments{
\item{data}{A dataframe of the attributes of microarray/RNA-seq data.}

\item{sample}{An optional numeric array that designates
the samples, that is rows, of the dataframe to be loaded.}
}
\value{
An object that represents the attributes of
    microarray/RNA-seq data.
    Internally, the object holds a dataframe whose number of entries
    (rows) equals that of the samples.
}
\description{
Converts attribute data from a dataframe into an object.
}

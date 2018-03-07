% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bGWAS_classFunctions.R
\name{extract_results_bGWAS}
\alias{extract_results_bGWAS}
\title{Exctract Results}
\usage{
extract_results_bGWAS(obj, SNPs = "significants")
}
\arguments{
\item{obj}{an object of class bGWAS}

\item{SNPs, }{"all" / "significants"}
}
\value{
a data.frame containing the results for all / significant SNPs
}
\description{
Extract results from an object of class bGWAS obtained
when using bGWAS() or bGWAS_fromPrior()
}
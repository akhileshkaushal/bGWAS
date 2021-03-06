% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bGWAS_classFunctions.R
\name{extract_results_bGWAS}
\alias{extract_results_bGWAS}
\title{Extract SNPs results from bGWAS results}
\usage{
extract_results_bGWAS(obj, SNPs = "significant")
}
\arguments{
\item{obj}{an object of class bGWAS}

\item{SNPs, }{"all" / "significant"}
}
\value{
a data.frame containing the results for all / significant SNPs
}
\description{
Extract SNPs results from an object of class bGWAS obtained
when using bGWAS() or bGWAS_fromPrior()
}

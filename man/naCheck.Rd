% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/naCheck.R
\name{naCheck}
\alias{naCheck}
\title{NA Checker}
\usage{
naCheck(physeq, trnl = TRUE, group = "Study")
}
\arguments{
\item{physeq}{phyloseq object of interest}

\item{trnl}{trnl or 12sv5}

\item{group}{grouping variable (e.g., study, location, etc.)}
}
\value{
for each study in the phyloseq object, returns the number of NA's present
}
\description{
Checks a phyloseq object for any taxa that have no superkingdom/kingdom assignment
}

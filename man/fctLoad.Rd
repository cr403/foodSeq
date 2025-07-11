% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fctLoad.R
\name{fctLoad}
\alias{fctLoad}
\title{Factor Loading Analysis}
\usage{
fctLoad(
  physeq,
  pca_output,
  PCx = 1,
  nTaxa = 10,
  name = "taxlabel",
  labWidth = 60,
  title = NULL,
  titleSize = 18,
  axisTitleSize = 16,
  textSize = 14
)
}
\arguments{
\item{physeq}{clr filtered phyloseq object on which pcaPlot() was run}

\item{pca_output}{output object from pcaPlot()}

\item{PCx}{principal component for x-axis of graph}

\item{nTaxa}{number of taxa to include on plot}

\item{name}{name of common name column}

\item{labWidth}{character length breakpoint for wrapping labels}

\item{title}{optional title}

\item{titleSize}{optional title size}

\item{axisTitleSize}{optional axis title size}

\item{textSize}{optional text size}
}
\value{
load.df = loadings data frame in descending order by PCx

plot = plot of top nTaxa taxa for PCx
}
\description{
This function generates factor loading plots.
}

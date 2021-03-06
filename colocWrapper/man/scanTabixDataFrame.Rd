% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/import_data.R
\name{scanTabixDataFrame}
\alias{scanTabixDataFrame}
\title{A general function to quickly import tabix indexed tab-separated files into data_frame}
\usage{
scanTabixDataFrame(tabix_file, param, ...)
}
\arguments{
\item{tabix_file}{Path to tabix-indexed text file}

\item{param}{A instance of GRanges, RangedData, or RangesList
provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()}

\item{...}{Additional parameters to be passed on to readr::read_delim()}
}
\value{
List of data_frames, one for each entry in the param GRanges object.
}
\description{
A general function to quickly import tabix indexed tab-separated files into data_frame
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CSCD.R
\name{CSCD}
\alias{CSCD}
\title{Performs gene expression decomposition.}
\usage{
CSCD(
  bulk.eset,
  sc.eset,
  min.p = NULL,
  markers = NULL,
  cell.types = "cellType",
  subj.names = "SubjectName",
  verbose = TRUE
)
}
\arguments{
\item{bulk.eset}{ExpressionSet with bulk data. Bulk RNA-seq data can be converted from a matrix with samples in columns and genes in rows to an ExpressionSet. See \code{\link{example_data}} for an example on how to create a bulk.eset object.}

\item{sc.eset}{ExpressionSet with single-cell data. Single-cell data requires additional information in the ExpressionSet, specifically cell-type labels and individual labels. See \code{\link{example_data}} for an example on how to create a sc.eset object.}

\item{min.p}{A percentage. This parameter is passed to the \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}} function (Butler et al., 2019) to pick the most relevant genes to each cell-type cluster. Users may pick a number between 0.3 and 0.5 for best results. The higher the value, the more genes to be excluded from the analysis.}

\item{markers}{A character vector containing marker genes
to be used in decomposition. If NULL provided, the method will use all available genes for decomposition.}

\item{cell.types}{Character string. A vector of cell-type labels.}

\item{subj.names}{Character string. A vector of individual labels that correspond to cells.}

\item{verbose}{Boolean. Whether to print log info during decomposition.
Errors will be printed regardless.}
}
\value{
A list. Slot \strong{bulk.props} contains a matrix of cell-type
proportion estimates with cell types as rows and individuals as columns.
Slot \strong{sc.props} contains a matrix of cell-type proportions
estimated directly from counting single-cell data. Slot \strong{transformed.bulk}
contains the covariance-based transformed bulk expression used for decomposition.
These values are generated by applying a linear transformation to the CPM expression.
Slot \strong{genes.used} contains a vector of genes used in decomposition. Slot
\strong{rnorm} contains Euclidean norm of the residuals for each individual's
proportion estimates.
}
\description{
Provides accurate cell-type proportion estimation by incorporating covariance structure
in given single-cell RNA-seq (scRNA-seq) and bulk RNA-seq datasets, see Karimnezhad (2022). The approach uses an extension
of the transformation used in Jew et al. (2020) implemented in the \code{\link[BisqueRNA:ReferenceBasedDecomposition]{BisqueRNA::ReferenceBasedDecomposition()}} function.
}
\references{
Butler, A. et al. (2019). Seurat: Tools for Single Cell Genomics. R package version, 4.1.1.

Jew, B. et al. (2020) Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nat Commun 11, 1971. https://doi.org/10.1038/s41467-020-15816-6

Jew, B. and Alvarez, M. (2020). BisqueRNA: Decomposition of Bulk Expression with Single-Cell Sequencing. R package version, 1.0.5.

Karimnezhad, A. (2022) More accurate estimation of cell composition in bulk expression through robust integration of single-cell information. https://doi.org/10.1101/2022.05.13.491858
}

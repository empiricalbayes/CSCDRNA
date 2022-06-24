#' computing the square root of a given high-dimension matrix using the Denman Beavers method
#'
#' Computes the square root of a given high-dimension matrix using the Denman Beavers method.
#' The original code is available here: https://github.com/XD-DENG/sqrt-matrix.
#' @examples
#' true_answer <- matrix(runif(100, 10, 100), 10, 10)
#' M <- true_answer %*% true_answer
#' inverse_matrix=DB_sqrtm(M)
#'
#' @param M a matrix
#' @return A matrix
#' @noRd
DB_sqrtm <- function(M){
  Y <- M
  Z <- diag(dim(M)[1])
  error <- 1
  error_limit <- 1.5e-8
  i=1
  while(error>error_limit){
    Y_old <- Y
    Y <- (Y_old+solve(Z))/2
    Z <- (Z+solve(Y_old))/2
    error <- max(abs(Y-Y_old))
    i=i+1
  }
  return(Y)
}

#' Find overlapping genes in single-cell data, bulk data, and marker genes.
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/reference_based.R
#'
#' @param sc.eset ExpressionSet with single-cell data
#' @param bulk.eset ExpressionSet with bulk data
#' @param markers Character vector. List of relevant marker genes
#' @param verbose Boolean. Print logging info
#' @return overlapping.genes Character vector. List of genes found in markers
#'   and both data sets.
#' @noRd
GetOverlappingGenes <- function(sc.eset, bulk.eset, markers, verbose) {
  bulk.genes <- Biobase::featureNames(bulk.eset)
  sc.genes <- Biobase::featureNames(sc.eset)
  overlapping.genes <- base::intersect(bulk.genes, sc.genes)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No overlapping genes found between bulk and ",
                            "single-cell expression."))
  }
  overlapping.genes <- base::intersect(overlapping.genes, markers)
  if (base::length(overlapping.genes) == 0) {
    base::stop(base::paste0("No marker genes found in both bulk and ",
                            "single-cell expression."))
  }
  if (verbose) {
    n.genes <- base::length(overlapping.genes)
    base::message(base::sprintf("Using %i genes in both", n.genes),
                  " bulk and single-cell expression.")
  }
  return(overlapping.genes)
}

#' Convert counts data in ExpressionSet to counts per million (CPM).
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/utils.R
#' @param eset ExpressionSet containing counts assay data.
#' @return eset ExpressionSet containing CPM assay data
#' @noRd
CountsToCPM <- function(eset) {
  Biobase::exprs(eset) <- base::sweep(Biobase::exprs(eset),
                                      2, base::colSums(Biobase::exprs(eset)),
                                      `/`) * 1000000
  indices <- base::apply(Biobase::exprs(eset), MARGIN=2,
                         FUN=function(column) {base::anyNA(column)})
  if (base::any(indices)) {
    n.cells <- base::sum(indices)
    base::stop(base::sprintf("Zero expression in selected genes for %i cells",
                             n.cells))
  }
  return(eset)
}

#' Remove genes in ExpressionSet with zero variance across samples.
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/utils.R
#'
#' @param eset ExpressionSet
#' @param verbose Boolean. Print logging info
#' @return eset ExpressionSet with zero variance genes removed
#' @noRd
FilterZeroVarianceGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, stats::var) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i zero variance genes.",
                                genes.filtered))
  }
  return(eset)
}

#' Remove genes in ExpressionSet with zero expression in all samples.
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/utils.R
#'
#' @param eset ExpressionSet
#' @param verbose Boolean. Print logging info
#' @return eset ExpressionSet with zero expression genes removed
#' @noRd
FilterUnexpressedGenes <- function(eset, verbose=TRUE) {
  indices <- (base::apply(Biobase::exprs(eset), 1, base::sum) != 0)
  indices <- indices & (! base::is.na(indices))
  if (base::sum(indices) < base::length(indices)) {
    eset <-
      Biobase::ExpressionSet(assayData=Biobase::exprs(eset)[indices,,drop=F],
                             phenoData=eset@phenoData)
  }
  if (verbose) {
    genes.filtered <- base::length(indices) - base::nrow(eset)
    base::message(base::sprintf("Filtered %i unexpressed genes.",
                                genes.filtered))
  }
  return(eset)
}

#' Generate reference profile for cell types identified in single-cell data.
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/reference_based.R
#'
#' Averages expression within each cell type across all samples to use as
#' reference profile.
#'
#' @param sc.eset ExpressionSet with single-cell data
#' @param cell.types A character string. Name of phenoData attribute in sc.eset
#'   that indicates cell type
#' @return sc.ref Matrix. Reference profile with number of gene rows by number
#'   of cell types columns.
#' @noRd
GenerateSCReference <- function(sc.eset, cell.types) {
  cell.labels <- base::factor(sc.eset[[cell.types]])
  all.cell.types <- base::levels(cell.labels)
  aggr.fn <- function(cell.type) {
    base::rowMeans(Biobase::exprs(sc.eset)[,cell.labels == cell.type, drop=F])
  }
  template <- base::numeric(base::nrow(sc.eset))
  sc.ref <- base::vapply(all.cell.types, aggr.fn, template)
  return(sc.ref)
}

#' Calculate cell proportions based on single-cell data.
#' The original code is available here: https://github.com/cran/BisqueRNA/blob/master/R/reference_based.R
#'
#' Returns proportion of each cell type out of total cells for each individual
#' in the single-cell ExpressionSet
#'
#' @param sc.eset ExpressionSet with single-cell data
#' @param subj.names A character string. Name of phenoData attribute in
#'   sc.eset that indicates individual ID.
#' @param cell.types A character string. Name of phenoData attribute in sc.eset
#'   that indicates cell type
#' @return sc.props Matrix. Cell proportions with number of cell types rows
#'   by number of individuals columns
#' @noRd
CalculateSCCellProportions <- function(sc.eset, subj.names, cell.types) {
  individual.labels <- base::factor(sc.eset[[subj.names]])
  individuals <- base::levels(individual.labels)
  cell.labels <- base::as.factor(sc.eset[[cell.types]])
  aggr.fn <- function(individual) {
    base::table(cell.labels[individual.labels == individual]) /
      base::length(cell.labels[individual.labels == individual])
  }
  sc.props <- base::sapply(individuals, aggr.fn)
  return(sc.props)
}


#' A function used internally in CSCD to remove genes appearing in more than one cluster.
#'
#' Returns a filtered matrix of the input matrix.
#'
#' @param M A matrix outputted by FindAllMarkers function of Seurat.
#' @noRd
FindAllMarkers_filter.by.cluster=function(M){
  a=NULL
  for (j in unique(M$gene)){
    a=rbind(a,c(j,length(M[M$gene==j,]$cluster)))
  }
  a=data.frame(a)
  colnames(a)<-c('gene','n.clusters')
  matrix.const=a[a$n.clusters==1,]
  M=M[matrix.const$gene,]
  return(M)
}


#' Performs gene expression decomposition.
#'
#' Provides accurate cell-type proportion estimation by incorporating covariance structure
#'in given single-cell RNA-seq (scRNA-seq) and bulk RNA-seq datasets, see Karimnezhad (2022). The approach uses an extension
#'of the transformation used in Jew et al. (2020) implemented in the \code{\link[BisqueRNA:ReferenceBasedDecomposition]{BisqueRNA::ReferenceBasedDecomposition()}} function.
#' @param bulk.eset ExpressionSet with bulk data. Bulk RNA-seq data can be converted from a matrix with samples in columns and genes in rows to an ExpressionSet. See \code{\link{example_data}} for an example on how to create a bulk.eset object.
#' @param sc.eset ExpressionSet with single-cell data. Single-cell data requires additional information in the ExpressionSet, specifically cell-type labels and individual labels. See \code{\link{example_data}} for an example on how to create a sc.eset object.
#' @param min.p A percentage. This parameter is passed to the \code{\link[Seurat:FindAllMarkers]{Seurat::FindAllMarkers()}} function (Butler et al., 2019) to pick the most relevant genes to each cell-type cluster. Users may pick a number between 0.3 and 0.5 for best results. The higher the value, the more genes to be excluded from the analysis.
#' @param markers A character vector containing marker genes
#'   to be used in decomposition. If NULL provided, the method will use all available genes for decomposition.
#' @param cell.types Character string. A vector of cell-type labels.
#' @param subj.names Character string. A vector of individual labels that correspond to cells.
#' @param verbose Boolean. Whether to print log info during decomposition.
#'   Errors will be printed regardless.
#' @return A list. Slot \strong{bulk.props} contains a matrix of cell-type
#'   proportion estimates with cell types as rows and individuals as columns.
#'   Slot \strong{sc.props} contains a matrix of cell-type proportions
#'   estimated directly from counting single-cell data. Slot \strong{transformed.bulk}
#'   contains the covariance-based transformed bulk expression used for decomposition.
#'   These values are generated by applying a linear transformation to the CPM expression.
#'   Slot \strong{genes.used} contains a vector of genes used in decomposition. Slot
#'   \strong{rnorm} contains Euclidean norm of the residuals for each individual's
#'   proportion estimates.
#' @importFrom methods new
#' @import Biobase
#' @import BisqueRNA
#' @import plyr
#' @import MAST
#' @import Seurat
#' @references Butler, A. et al. (2019). Seurat: Tools for Single Cell Genomics. R package version, 4.1.1.
#' @references Jew, B. et al. (2020) Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nat Commun 11, 1971. https://doi.org/10.1038/s41467-020-15816-6
#' @references Jew, B. and Alvarez, M. (2020). BisqueRNA: Decomposition of Bulk Expression with Single-Cell Sequencing. R package version, 1.0.5.
#' @references Karimnezhad, A. (2022) More accurate estimation of cell composition in bulk expression through robust integration of single-cell information. https://doi.org/10.1101/2022.05.13.491858
#'
#' @export
CSCD<-function(bulk.eset, sc.eset,min.p=NULL,markers=NULL,
               cell.types="cellType",subj.names="SubjectName",verbose=TRUE) {
  old.cpm=FALSE
  if ((! methods::is(sc.eset, "ExpressionSet")) ||
      (! methods::is(bulk.eset, "ExpressionSet"))) {
    base::stop("Expression data should be in ExpressionSet")
  }
  else if (! cell.types %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Cell-type label \"%s\" ", cell.types),
               "not found in single-cell ExpressionSet varLabels.")
  }
  else if (! subj.names %in% Biobase::varLabels(sc.eset)) {
    base::stop(base::sprintf("Individual label \"%s\"", subj.names),
               " not found in single-cell ExpressionSet varLabels.")
  }
  n.sc.individuals <-
    base::length(base::levels(base::factor(sc.eset[[subj.names]])))
  if (n.sc.individuals == 1) {
    base::stop("Only one individual detected in single-cell data. At least ",
               "two subjects are needed (three or more recommended).")
  }
  else if (n.sc.individuals == 2) {
    base::warning("Only two individuals detected in single-cell data. While ",
                  "CSCD will run, we recommend at least three subjects for",
                  " reliable performance.")
  }
  n.cell.types <-
    base::length(base::levels(base::factor(sc.eset[[cell.types]])))
  if (n.cell.types == 1) {
    base::stop("Single-cell pheno data indicates only one cell type",
               " present. No need for decomposition.")
  }
  if (verbose) {
    base::message(base::sprintf("Decomposing into %i cell types.",
                                n.cell.types))
  }
  if (base::is.null(markers)) {
    markers <- Biobase::featureNames(sc.eset)
  }
  else {
    markers <- base::unique(base::unlist(markers))
  }


  if (base::is.null(min.p)) {
    min.p <- 0.5
  }

  test.use =  'MAST'# 'bimod'#"roc"#
  sc_raw=exprs(sc.eset)
  sc <- Seurat::CreateSeuratObject(counts = sc_raw)
  sc$replicate <- sc.eset$SubjectName
  Idents(sc)<-as.character(sc.eset$cellType)
  id=unique(sc.eset$cellType)

  markers.identified=Seurat::FindAllMarkers(object = sc,test.use = test.use,only.pos = T,min.diff.pct = min.p,min.pct = 0,pseudocount.use = .01)
  markers.identified=FindAllMarkers_filter.by.cluster(markers.identified)
  g.w.PT=rownames(markers.identified)
  sc.eset=sc.eset[which(rownames(sc.eset)%in%g.w.PT),]

  genes <- GetOverlappingGenes(sc.eset, bulk.eset, markers, verbose)
  if (verbose) {
    base::message("Converting single-cell counts to CPM and ",
                  "filtering zero variance genes.")
  }


  sc.eset <- CountsToCPM(sc.eset)
  sc.eset <-Biobase::ExpressionSet(assayData=Biobase::exprs(sc.eset)[genes,],
                                   phenoData=sc.eset@phenoData)
  sc.eset <- FilterZeroVarianceGenes(sc.eset, verbose)
  if (verbose) {
    base::message("Converting bulk counts to CPM and filtering",
                  " unexpressed genes.")
  }
  bulk.eset <- CountsToCPM(bulk.eset)
  bulk.eset <- Biobase::ExpressionSet(assayData=Biobase::exprs(bulk.eset)[genes,],
                                      phenoData=bulk.eset@phenoData)

  bulk.eset <- FilterUnexpressedGenes(bulk.eset, verbose)
  genes <- base::intersect(Biobase::featureNames(sc.eset),
                           Biobase::featureNames(bulk.eset))
  if (base::length(genes) == 0) {
    base::stop("Zero genes remaining after filtering and ",
               "intersecting bulk, single-cell, and marker genes.")
  }
  if (verbose) {
    n.cells <- base::ncol(sc.eset)
    base::message("Generating single-cell based reference from ",
                  sprintf("%i cells.\n", n.cells))
  }
  sc.ref <- GenerateSCReference(sc.eset, cell.types)[genes,,drop=F]
  sc.props <- CalculateSCCellProportions(sc.eset, subj.names, cell.types)
  sc.props <- sc.props[base::colnames(sc.ref),,drop=F]

  Z <- GenerateSCReference(sc.eset, cell.types)[genes,,drop=F]
  p <- CalculateSCCellProportions(sc.eset, subj.names, cell.types)
  p <- p[base::colnames(Z),,drop=F]

  Y <- Z %*% p
  X <- Biobase::exprs(bulk.eset)[genes,,drop=F]
  sample.names <- base::colnames(Biobase::exprs(bulk.eset))
  template <- base::numeric(base::length(sample.names))
  base::names(template) <- sample.names
  base::message('applying the covariance-based transformation ...')
  cov.Y.multiv.bisque<-nlshrink::linshrink_cov(t(Y))
  cov.Y.sqrt.multiv.bisque<-DB_sqrtm(cov.Y.multiv.bisque)
  rownames(cov.Y.multiv.bisque) <-rownames(cov.Y.sqrt.multiv.bisque) <- rownames(X)
  colnames(cov.Y.multiv.bisque) <-colnames(cov.Y.sqrt.multiv.bisque) <- rownames(X)
  cov.X<-nlshrink::linshrink_cov(t(X))
  cov.X.inv<-solve(cov.X)
  cov.X.inv.sqrt<-DB_sqrtm(cov.X.inv)
  rownames(cov.X.inv) <-rownames(cov.X.inv.sqrt) <- rownames(X)
  colnames(cov.X.inv) <-colnames(cov.X.inv.sqrt) <- rownames(X)
  X.star.multiv<-((cov.Y.sqrt.multiv.bisque%*%cov.X.inv.sqrt)%*%(X-rowMeans(X)))*sqrt((n.sc.individuals-1) /(n.sc.individuals+1)) + rowMeans(Y)
  base::rownames(X.star.multiv) <- base::rownames(Z)
  base::colnames(X.star.multiv) <- sample.names
  X.star.multiv[X.star.multiv<0]=0
  base::message('estimating cell-type proportions ...')
  E <- base::matrix(1, nrow = n.cell.types, ncol = n.cell.types)
  f <- base::rep(1, n.cell.types)
  G <- base::diag(n.cell.types)
  h <- base::rep(0, n.cell.types)
  res.multiv.bisque <- base::as.matrix(base::apply(t( X.star.multiv), 1,
                                                   function(b) {
                                                     sol <- limSolve::lsei(Z, b,
                                                                           E, f, G, h)
                                                     sol.p <- sol$X
                                                     sol.r <- base::sqrt(sol$solutionNorm)
                                                     return(base::append(sol.p, sol.r))
                                                   }))

  rownames( res.multiv.bisque )[n.cell.types+1]<-'rnorm'
  results=list(bulk.props=res.multiv.bisque[1:n.cell.types,],sc.props=p,
               transformed.bulk=X.star.multiv,rnorm=res.multiv.bisque['rnorm',],genes.used=genes)
  return(results)
}


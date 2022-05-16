# CSCDRNA
 Covariance Based Single-Cell Decomposition of Bulk Expression Data
 
 ### Cell-type decomposition
This approach provides accurate cell-type proportion estimation by incorporating covariance structure in a given set of single-cell RNA-seq (scRNA-seq) and bulk RNA-seq datasets, see Karimnezhad (2022). The approach uses an extension of the transformation used in Jew et al. (2020).

## Installation

The package can be installed from the GitHub repository
```r
devtools::install_github("empiricalbayes/CSCDRNA")
```

## Getting Started
You can load Bisque as follows:

```r
library(CSCDRNA)
```

## How to run CSCD on an example data containing synthetic bulk and single-cell datasets.
```r
#Load example data.
data(example_data)

#Build ExpressionSet with bulk data.
bulk.eset <- Biobase::ExpressionSet(assayData = example_data$bulk.matrix)

#Build ExpressionSet with single-cell data.
sc.counts.matrix=example_data$sc.counts.matrix
individual.labels=example_data$individual.labels
cell.type.labels=example_data$cell.type.labels
sample.ids <- colnames(sc.counts.matrix)
#individual.labels and cell.types should be in the same order as in sample.ids.

sc.pheno <- data.frame(check.names=FALSE, check.rows=FALSE,
                        stringsAsFactors=FALSE,row.names=sample.ids,
                        SubjectName=individual.labels,cellType=cell.type.labels)
 sc.meta <- data.frame(labelDescription=c("SubjectName","cellType"),
                       row.names=c("SubjectName","cellType"))
sc.pdata <- new("AnnotatedDataFrame",data=sc.pheno, varMetadata=sc.meta)
sc.eset <- Biobase::ExpressionSet(assayData=sc.counts.matrix,phenoData=sc.pdata)

#Run CSCD on the example data.
analysis <- CSCD(bulk.eset=bulk.eset,sc.eset= sc.eset,
                 min.p=0.3,markers=NULL,cell.types="cellType",
                  subj.names="SubjectName",verbose=TRUE)


#Estimated cell-type proportions.
analysis$bulk.props

#Cell-type proportions estimated directly by counting single-cell data.
analysis$sc.props

#The covariance based transformed bulk expression used for decomposition.
analysis$transformed.bulk.

#Genes used in the decomposition.
analysis$genes.used

#Euclidean norm of the residuals for each individual's proportion estimates.
analysis$rnorm
```


## References
Jew, B. et al. (2020) Accurate estimation of cell composition in bulk expression through robust integration of single-cell information. Nat Commun 11, 1971. https://doi.org/10.1038/s41467-020-15816-6

Karimnezhad, A. (2022) More accurate estimation of cell composition in bulk expression through robust integration of single-cell information. https://www.biorxiv.org/content/10.1101/2022.05.13.491858v1

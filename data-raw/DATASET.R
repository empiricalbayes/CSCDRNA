## code to prepare `DATASET` dataset goes here
sim.data <- BisqueRNA::SimulateData(n.ind=5, n.genes=500, n.cells=4,
                                    cell.types=c("CT1","CT2","CT3","CT4"),
                                    avg.props= c(0.4,0.3,0.2,0.1))
example_data=NULL
example_data$bulk.matrix=Biobase::exprs(sim.data$bulk.eset)
example_data$sc.counts.matrix=Biobase::exprs(sim.data$sc.eset)
example_data$individual.labels=sim.data$sc.eset$SubjectName
example_data$cell.type.labels=sim.data$sc.eset$cellType
example_data$sample.ids <- colnames(example_data$sc.counts.matrix)

usethis::use_data(example_data, overwrite = TRUE)

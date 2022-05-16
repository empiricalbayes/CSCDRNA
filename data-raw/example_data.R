## code to prepare `example_data` dataset goes here
set.seed(100)
sim.data <- BisqueRNA::SimulateData(n.ind=5, n.genes=100, n.cells=4,
                                    cell.types=c("Endo","Podo","PT","DCT"),
                                    avg.props= c(0.55,0.25,0.15,0.05))
example_data=NULL
example_data$bulk.matrix=Biobase::exprs(sim.data$bulk.eset)
example_data$sc.counts.matrix=Biobase::exprs(sim.data$sc.eset)
example_data$individual.labels=sim.data$sc.eset$SubjectName
example_data$cell.type.labels=sim.data$sc.eset$cellType
example_data$sample.ids <- colnames(example_data$sc.counts.matrix)

usethis::use_data(example_data, overwrite = TRUE)

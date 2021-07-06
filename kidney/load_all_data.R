setwd("~/Documents/main_files/Explain/gencoder/project_gencoder/")

kidney_full <- zellkonverter::readH5AD("./data/external/kidney/Mature_Full_v3.h5ad")

count_kidney_full <- kidney_full@assays@data$counts
meta_kidney_full.main <- do.call('cbind',kidney_full@colData@listData)

meta_kidney_full <- data.frame("Short_Sample" = as.factor(meta_kidney_full.main[,2]),
                               "Project" = as.factor(meta_kidney_full.main[,3]),
                               "Experiment" = as.factor(meta_kidney_full.main[,4]),
                               "celltype" = as.factor(meta_kidney_full.main[,5]),
                               "compartment" = as.factor(meta_kidney_full.main[,6]),
                               "broad_celltype" = as.factor(meta_kidney_full.main[,7]))

meta_kidney_full <- mltools::one_hot(data.table::data.table(meta_kidney_full))

save(list = c("count_kidney_full","meta_kidney_full","meta_kidney_full.main","HVG_parenchyma"),file = "./data/internal/kidney/Mature_Full_v3.Rdata")



kidney_parenchyma <- zellkonverter::readH5AD("./data/external/kidney/Mature_non_PT_parenchyma_v2.1.h5ad")

count_kidney_parenchyma <- kidney_parenchyma@assays@data$X

HVG_parenchyma <- row.names(count_kidney_parenchyma)[order(apply(count_kidney_parenchyma,1,var),decreasing = T)]

meta_kidney_parenchyma.main <- do.call('cbind',kidney_parenchyma@colData@listData)


meta_kidney_parenchyma <- data.frame("celltype.parenchyma" = as.factor(meta_kidney_parenchyma.main[,1]),
                               "broad_celltype.parenchyma" = as.factor(meta_kidney_parenchyma.main[,2]),
                               "Experiment.parenchyma" = as.factor(meta_kidney_parenchyma.main[,3]),
                               "Project.parenchyma" = as.factor(meta_kidney_parenchyma.main[,4]))

meta_kidney_parenchyma <- mltools::one_hot(data.table::data.table(meta_kidney_parenchyma))

save(list = c("count_kidney_parenchyma","meta_kidney_parenchyma.main","meta_kidney_parenchyma","HVG_parenchyma"),file = "./data/internal/kidney/Mature_non_PT_parenchyma_v2.1.Rdata")

write.csv(levels(kidney_parenchyma@colData@listData[["celltype"]]),"./data/workflow/main_kidney_list.csv",quote=F,row.names=T,col.names=1)


kidney_immune <- zellkonverter::readH5AD("./data/external/kidney/Mature_Immune_v2.1.h5ad")

count_kidney_immune <- kidney_immune@assays@data$X

meta_kidney_immune.main <- do.call('cbind',kidney_immune@colData@listData)


meta_kidney_immune <- data.frame("celltype.immune" = as.factor(meta_kidney_immune.main[,1]),
                                     "broad_celltype.immune" = as.factor(meta_kidney_immune.main[,2]),
                                     "Experiment.immune" = as.factor(meta_kidney_immune.main[,3]),
                                     "Project.immune" = as.factor(meta_kidney_immune.main[,4]))

meta_kidney_immune <- mltools::one_hot(data.table::data.table(meta_kidney_immune))

save(list = c("count_kidney_immune","meta_kidney_immune.main","meta_kidney_immune"),file = "./data/internal/kidney/Mature_Immune_v2.1.Rdata")






fname <- "./data/external/kidney/GSE109564_Kidney.biopsy.dge.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE109564&format=file&file=GSE109564%5FKidney%2Ebiopsy%2Edge%2Etxt%2Egz",fname) }
GSE109564_Kidney <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE109564_Kidney"),file = "./data/internal/kidney/GSE109564_Kidney.Rdata")


fname <- "./data/external/kidney/GSE160048_human_glom_single_cell_rpkms.txt.gz"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE160048&format=file&file=GSE160048%5Fhuman%5Fglom%5Fsingle%5Fcell%5Frpkms%2Etxt%2Egz",fname) }
GSE160048_human_glom <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE160048_human_glom"),file = "./data/internal/kidney/GSE160048_human_glom.Rdata")


fname <- "./data/external/kidney/GSE131685_RAW.tar"
if (!file.exists(fname)) { download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131685&format=file",fname) }
GSE131685_RAW <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE131685_RAW"),file = "./data/internal/kidney/GSE131685_RAW.Rdata")

# biopsy.seurat.Robj is at:
## https://data.humancellatlas.org/explore/projects/027c51c6-0719-469f-a7f5-640fe57cbece/expression-matrices
load("./data/external/kidney/biopsy.seurat.Robj",verbose = T)
biopsy <- Seurat::UpdateSeuratObject(biopsy)
kidney.counts <- biopsy@assays$RNA@counts
colnames(kidney.counts) <- biopsy@active.ident
#save(list = c("kidney.counts"),file = "./data/internal/kidney/kidney.counts.Rdata")



fname <- "./data/external/kidney/GSE118184_BDNF_inhibitor.dge.txt.gz"
GSE118184_BDNF_inhibitor <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_BDNF_inhibitor"),file = "./data/internal/kidney/GSE118184_BDNF_inhibitor.Rdata")


fname <- "./data/external/kidney/GSE118184_Human_kidney_snRNA.dge.txt.gz"
GSE118184_Human_kidney_snRNA <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Human_kidney_snRNA"),file = "./data/internal/kidney/GSE118184_Human_kidney_snRNA.Rdata")


fname <- "./data/external/kidney/GSE118184_Morizane_ES.dge.txt.gz"
GSE118184_Morizane_ES <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Morizane_ES"),file = "./data/internal/kidney/GSE118184_Morizane_ES.Rdata")


fname <- "./data/external/kidney/GSE118184_Morizane_iPS_Batch1_2.dge.txt.gz"
GSE118184_Morizane_iPS_Batch1_2 <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Morizane_iPS_Batch1_2"),file = "./data/internal/kidney/GSE118184_Morizane_iPS_Batch1_2.Rdata")


fname <- "./data/external/kidney/GSE118184_Takasato_ES.dge.txt.gz"
GSE118184_Takasato_ES <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Takasato_ES"),file = "./data/internal/kidney/GSE118184_Takasato_ES.Rdata")


fname <- "./data/external/kidney/GSE118184_Takasato_iPS_Batch1_2.dge.txt.gz"
GSE118184_Takasato_iPS_Batch1_2 <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Takasato_iPS_Batch1_2"),file = "./data/internal/kidney/GSE118184_Takasato_iPS_Batch1_2.Rdata")


fname <- "./data/external/kidney/GSE118184_Takasato_iPS_Batch3.dge.txt.gz"
GSE118184_Takasato_iPS_Batch3 <- read.table(fname, header=TRUE, row.names=1, check.names=FALSE)
#save(list = c("GSE118184_Takasato_iPS_Batch3"),file = "./data/internal/kidney/GSE118184_Takasato_iPS_Batch3.Rdata")


spatial_nanostring <- read.delim("./data/external/kidney/Kidney_Raw_TargetCountMatrix.txt",row.names = 1)
spatial_nanostring <- spatial_nanostring[,order(colnames(spatial_nanostring),decreasing = F)]
pixel_nanostring <- read.csv("./data/external/kidney/TRAIN_pixels_color_32.csv",header=T,row.names=1)
unrelated_pixel_nanostring <- read.csv("./data/external/kidney/TEST_pixels_color_32.csv",header=T,row.names=1)

save(list = c("pixel_nanostring","unrelated_pixel_nanostring","spatial_nanostring"),file = "./data/internal/kidney/spatial_nanostring.Rdata")

HVG_spatial_nanostring <- row.names(spatial_nanostring)[order(apply(spatial_nanostring,1,var),decreasing = T)]

main_HVG <- Reduce('intersect',list(
  HVG_parenchyma,
  HVG_spatial_nanostring
))

side_HVG <- Reduce('intersect',list(
  HVG_parenchyma,
  HVG_spatial_nanostring,
  row.names(GSE118184_BDNF_inhibitor),
  row.names(GSE118184_Human_kidney_snRNA),
  row.names(GSE118184_Morizane_ES),
  row.names(GSE118184_Morizane_iPS_Batch1_2),
  row.names(GSE118184_Takasato_ES),
  row.names(GSE118184_Takasato_iPS_Batch1_2),
  row.names(GSE118184_Takasato_iPS_Batch3)
))


save(list=c("HVG_spatial_nanostring","HVG_parenchyma","main_HVG","side_HVG"),file = "./data/internal/kidney/HVG.Rdata")



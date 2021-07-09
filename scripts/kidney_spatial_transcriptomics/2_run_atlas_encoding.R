load("./data/internal/kidney/Mature_Immune_v2.1.Rdata")
load("./data/internal/kidney/Mature_non_PT_parenchyma_v2.1.Rdata")
load(file = "./data/internal/kidney/Mature_fetal_v3.Rdata")
load("./data/internal/kidney/HVG.Rdata")

devtools::load_all("~/Documents/main_files/GCP/explain-dpm/clean/gcproc_generalised_canonical_procrustes_2020_06_04/gcproc/R/")

n_genes = length(main_HVG)

load("./data/internal/kidney/spatial_nanostring.Rdata")

y <- Matrix::t(Matrix::Matrix(as.matrix(spatial_nanostring[main_HVG[1:n_genes],]),sparse = T))
x <- Matrix::t(count_kidney_parenchyma[main_HVG[1:n_genes],])
z <- Matrix::t(count_kidney_fetal[main_HVG[1:n_genes],])
p <- pixel_nanostring
r <- unrelated_pixel_nanostring

fetal_atlas.gcproc <- gcproc::gcproc(x=z,
                               y=z,
                               config = list(k_dim = 70,
                                             j_dim = 70,
                                             eta=1e-2,
                                             max_iter=10000,
                                             min_iter = 15,
                                             tol=1e-2,
                                             log=F,
                                             center=T,
                                             scale.z=T,
                                             batches=128,
                                             cores=8,
                                             verbose=T,
                                             init="svd"),
                               seed = 1)


save(fetal_atlas.gcproc,"./data/workflow/fetal_atlas.Rdata")
rm(fetal_atlas.gcproc)


atlas.gcproc <- gcproc::gcproc(x=x,
                               y=x,
                               config = list(k_dim = 70,
                                                     j_dim = 70,
                                                     eta=0.01,
                                                     max_iter=1500,
                                                     min_iter = 15,
                                                     tol=1e-5,
                                                     log=F,
                                                     center=T,
                                                     scale.z=T,
                                                     batches=512,
                                                     cores=8,
                                                     verbose=T,
                                                     init="svd"),
                               seed = 1)




spatial_nanostring.gcproc <- gcproc::gcproc(x=x,
                               y=y,
                               config = list(k_dim = 70,
                                             j_dim = 70,
                                             eta=0.01,
                                             max_iter=1500,
                                             min_iter = 15,
                                             tol=1e-5,
                                             log=F,
                                             center=T,
                                             scale.z=T,
                                             batches=64,
                                             cores=8,
                                             verbose=T,
                                             init="svd"),
                               seed = 1,
                               anchor = list(
                                 anchor_y.sample = NULL,
                                 anchor_y.feature = NULL,
                                 anchor_x.sample = atlas.gcproc$main.parameters$alpha.L.J,
                                 anchor_x.feature = NULL
                               )
                               )



pixel_nanostring.gcproc <- gcproc::gcproc(x=p,
                                            y=y,
                                            config = list(k_dim = 70,
                                                          j_dim = 70,
                                                          eta=0.01,
                                                          max_iter=1500,
                                                          min_iter = 15,
                                                          tol=1e-8,
                                                          log=F,
                                                          center=T,
                                                          scale.z=T,
                                                          batches=64,
                                                          cores=8,
                                                          verbose=T,
                                                          init="svd"),
                                            seed = 1,
                                          anchor = list(
                                            anchor_y.sample = spatial_nanostring.gcproc$main.parameters$alpha.L.K,
                                            anchor_y.feature = NULL,
                                            anchor_x.sample = NULL,
                                            anchor_x.feature = NULL
                                          )
)




unrelated_pixel_nanostring.gcproc <- gcproc::gcproc(x=r,
                                            y=x,
                                            config = list(k_dim = 70,
                                                          j_dim = 70,
                                                          eta=0.01,
                                                          max_iter=1500,
                                                          min_iter = 15,
                                                          tol=1e-8,
                                                          log=F,
                                                          center=T,
                                                          scale.z=T,
                                                          batches=64,
                                                          cores=8,
                                                          verbose=T,
                                                          init="svd"),
                                            seed = 1,
                                            anchor = list(
                                              anchor_y.sample = NULL,
                                              anchor_y.feature = spatial_nanostring.gcproc$main.parameters$u.beta,
                                              anchor_x.sample = NULL,
                                              anchor_x.feature = NULL
                                            )
                                          )


v <- unrelated_pixel_nanostring.gcproc$main.parameters$v.beta[-c(1:dim(meta_kidney)[1]),]
u <- unrelated_pixel_nanostring.gcproc$main.parameters$u.beta
  
u.inv <- MASS::ginv(t(u))
v.inv <- MASS::ginv(v%*%t(v))

gene2pixel <- v%*%t(u)%*%u.inv
pixel2gene <- u%*%t(v)%*%v.inv

colnames(pixel2gene) <- row.names(gene2pixel) <- c(main_HVG[1:n_genes])
colnames(x)[-c(1:dim(meta_kidney)[1])] <- c(main_HVG[1:n_genes])



write.csv(round(gene2pixel,2),"./data/workflow/v6/gene2pixel.csv",quote=F,row.names=T,col.names=T)
write.csv(round(pixel2gene,2),"./data/workflow/v6/pixel2gene.csv",quote=F,row.names=T,col.names=T)

write.csv(round(u,5),"./data/workflow/v6/pixel_reduction.csv",quote=F,row.names=T,col.names=T)
write.csv(round(v,5),"./data/workflow/v6/gene_reduction.csv",quote=F,row.names=T,col.names=T)

write.csv(round(as.matrix(x),2),"./data/workflow/v6/main_image_as_gene_v1.csv",quote=F,row.names=T,col.names=T)

save.image(file = "./data/workflow/v6/all.nanostring_via_atlas.gcproc.Rdata")


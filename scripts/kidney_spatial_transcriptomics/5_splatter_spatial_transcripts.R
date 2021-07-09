load("./data/workflow/fetal_atlas.Rdata")

set.seed(1)
cell_fetal_id <- do.call('c',lapply(unique(meta_kidney_fetal.main[,8]),function(X){
  sample(which(meta_kidney_fetal.main[,8]%in%X),size = 100)
}))

cell_fetal_cell_type <- rep(fetal_cell_type_labels[unique(meta_kidney_fetal.main[,8])],each=100)


library(splatter)
x <- count_kidney_fetal[main_HVG[1:1000],cell_fetal_id]
params <- splatter::splatEstimate(as.matrix(x))
sim <- splatter::splatSimulate(params, nGenes = 1000)

y <- gcproc::prepare_data(sim@assays@data$counts,rnorm(10),log = F,center = T,scale.z = T)
y <- y$x


n_gene = 10000

x <- count_kidney_fetal[main_HVG[1:n_gene],cell_fetal_id]
K <- atlas.gcproc$main.parameters$alpha.L.K[,cell_fetal_id]
v <- atlas.gcproc$main.parameters$v.beta[c(1:n_gene),]
J <- atlas.gcproc$main.parameters$alpha.L.J[,cell_fetal_id]


x.rd <- K%*%t(x)%*%v
y.rd <- J%*%t(y)
u <- as.matrix(MASS::ginv(t(y.rd)%*%y.rd)%*%t(y.rd)%*%x.rd)

new.y.rd <- J%*%t(y)%*%u
beta <- as.matrix(MASS::ginv(t(new.y.rd)%*%new.y.rd)%*%t(new.y.rd)%*%K%*%t(x))

new.y <- t(t(y)%*%u%*%beta)
row.names(new.y) <- main_HVG[1:10000]




library(ggplot2)
library(grid)
library(gridExtra)


for (k in c(1:50)){
  index <- k
  main_new.y <- t(new.y)[1:4,]
  
  if (k < 25){
    main_new.y[,"WT1"] <- main_new.y[,"WT1"]*k+0.3*k
    main_new.y[,"PODXL"] <- main_new.y[,"PODXL"]*k+0.3*k
    main_new.y[,"NPHS1"] <- main_new.y[,"NPHS1"]*k+0.3*k
    main_new.y[,"NPHS2"] <- main_new.y[,"NPHS2"]*k+0.3*k
  }
  if (k >= 25){
    k <- k - 24
    main_new.y[,"VEGFA"] <- main_new.y[,"VEGFA"]*k+0.3*k
    main_new.y[,"VIM"] <- main_new.y[,"VIM"]+k+0.3*k
    main_new.y[,"CDH1"] <- main_new.y[,"CDH1"]+k+0.3*k
    main_new.y[,"AQP2"] <- main_new.y[,"AQP2"]+k+0.3*k
    
    # Provide continuity from discontinuous if statement
    k <- 25
    main_new.y[,"WT1"] <- main_new.y[,"WT1"]*k+0.3*k
    main_new.y[,"PODXL"] <- main_new.y[,"PODXL"]*k+0.3*k
    main_new.y[,"NPHS1"] <- main_new.y[,"NPHS1"]*k+0.3*k
    main_new.y[,"NPHS2"] <- main_new.y[,"NPHS2"]*k+0.3*k
    
  }

  new.tiles.y <- main_new.y%*%pixel2gene[1:n_gene,]
  
  
  convert_to_RGB <- function(X){
    x <- ecdf(X)
    j <- x(X)
    X <- array(j,dim=dim(X))
    return(X)
  }
  
  # for (i in 1:4){
  #   final_cell_image <- convert_to_RGB(aperm(array((new.tiles.y[i,]),dim=c(3,32,32)),c(2,3,1)))
  #   g <- rasterGrob(final_cell_image, interpolate=TRUE)
  #   a <- qplot(1:10, 1:10, geom="blank") +
  #     annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
  #   plot(a)
  # }
  
  
  final_cell_image1 <- convert_to_RGB(aperm(array((new.tiles.y[1,]),dim=c(3,32,32)),c(2,3,1)))
  final_cell_image2 <- convert_to_RGB(aperm(array((new.tiles.y[2,]),dim=c(3,32,32)),c(2,3,1)))
  final_cell_image3 <- convert_to_RGB(aperm(array((new.tiles.y[3,]),dim=c(3,32,32)),c(2,3,1)))
  final_cell_image4 <- convert_to_RGB(aperm(array((new.tiles.y[4,]),dim=c(3,32,32)),c(2,3,1)))
  
  final_cell_image1.2 <- array(0,dim=c(32,48,3))
  final_cell_image1.2[,1:16,] <- final_cell_image1[,1:16,]
  final_cell_image1.2[,32:48,] <- final_cell_image2[,16:32,]
  
  for (i in 1:16){
    alpha <- (16-i)/16
    final_cell_image1.2[,16+i,] <- (final_cell_image1[,16+i,]*alpha+(1-alpha)*final_cell_image2[,i,])
  }
  
  
  
  final_cell_image3.4 <- array(0,dim=c(32,48,3))
  final_cell_image3.4[,1:16,] <- final_cell_image3[,1:16,]
  final_cell_image3.4[,32:48,] <- final_cell_image4[,16:32,]
  
  for (i in 1:16){
    alpha <- (16-i)/16
    final_cell_image3.4[,16+i,] <- (final_cell_image3[,16+i,]*alpha+(1-alpha)*final_cell_image4[,i,])
  }
  
  
  
  
  
  final_cell_image1.2.3.4 <- array(0,dim=c(48,48,3))
  final_cell_image1.2.3.4[1:16,,] <- final_cell_image1.2[1:16,,]
  final_cell_image1.2.3.4[32:48,,] <- final_cell_image3.4[16:32,,]
  
  for (i in 1:16){
    alpha <- (16-i)/16
    final_cell_image1.2.3.4[16+i,,] <- (final_cell_image1.2[16+i,,]*alpha+(1-alpha)*final_cell_image3.4[i,,])
  }
  
  
  g <- rasterGrob(final_cell_image1.2.3.4, interpolate=TRUE)
  a <- qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

  
  gene_vector_i.j <- c()
  for (i in 1:16){
    for (j in 1:16){
      gene_vector_i.j <- cbind(gene_vector_i.j,(c(final_cell_image1.2.3.4[i:(i+31),j:(j+31),])))
    }
  }
  
  
  gene_grid <- t(gene_vector_i.j)%*%gene2pixel
  colnames(gene_grid) <- main_HVG
  
  
  gs <- list(a)
  gs <- c(gs,lapply(c("WT1","PODXL","SIX2","CDH1","AQP2","VIM"),function(X){
    
    g <- rasterGrob(convert_to_RGB(array(gene_grid[,X],dim=c(16,16))), interpolate=TRUE)
    a <- qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+ggtitle(X)
    
  }))
  
  lay <- rbind(c(4,1,1,1,5),
               c(2,1,1,1,6),
               c(3,1,1,1,7))
  
  final.p <- gridExtra::arrangeGrob(grobs = gs, layout_matrix = lay)
  
  ggsave(plot = final.p,device = "png",filename = paste("./figures/splatter_visualised/main___",index,".png",sep=""),width = 15,height = 9)
  
}

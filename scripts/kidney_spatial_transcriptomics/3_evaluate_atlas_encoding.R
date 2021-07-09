
convert_to_RGB <- function(x){
  x[x<0] <- 0
  x[x>255] <- 255
  x <- x/255
  return(x)
}


convert_to_HEAT <- function(x){
  e <- ecdf(x)
  j <- e(x)
  x <- array(j,dim=dim(x))
  return(x)
}

load(file = "./data/workflow backup/all.nanostring_via_atlas.gcproc.Rdata")

main_kidney_list <- read.csv("./data/workflow/v5/main_kidney_list.csv")
load("./data/internal/kidney/Mature_non_PT_parenchyma_v2.1.Rdata")
load("./data/internal/kidney/Mature_Immune_v2.1.Rdata")


all_list <- c(meta_kidney_parenchyma.main[,1],meta_kidney_immune.main[,1])

gene_data <- x[,-c(1:dim(meta_kidney)[1])]
top_gene_var_ID <- order(apply(gene_data,2,var),decreasing=T)
top_IDS_pod <- intersect(top_gene_var_ID,which(all_list==1))
top_IDS_tub <- intersect(top_gene_var_ID,which(all_list==4))

gene_tiles <- t((unrelated_pixel_nanostring.gcproc[["transformed.data"]][["x"]])%*%(pixel2gene))
row.names(gene_tiles) <- colnames(gene_data)



g <- (convert_to_HEAT(aperm(array((pixel_nanostring.gcproc$main.parameters$u.beta[,1]),dim=c(3,32,32)),c(2,3,1))))
g <- rasterGrob(g, interpolate=TRUE)
p13 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
p13


g <- (convert_to_HEAT(aperm(array((unrelated_pixel_nanostring.gcproc$main.parameters$u.beta[,1]),dim=c(3,32,32)),c(2,3,1))))
g <- rasterGrob(g, interpolate=TRUE)
p13 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
  annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
p13



if (F){
  library(png)
  library(ggplot2)
  library(grid)
  
  img <- png::readPNG("./data/workflow/v5/normal4_eg_1280.png")
  g <- rasterGrob(img, interpolate=TRUE)
  
  
  p1 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["SIX1",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p2 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("SIX1")
  p2
  
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["SIX2",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p3 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("SIX2")
  p3
  
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["HLA-A",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p4 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("HLA-A")
  p4
  
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["NPHS1",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p5 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("NPHS1")
  p5
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["NPHS2",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p6 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("NPHS2")
  p6
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["HLA-B",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p7 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("HLA-B")
  p7
  
  
  library(reshape2)
  g <- melt(convert_to_HEAT(t(apply(matrix(((gene_tiles["HLA-C",])),nrow=40,ncol=40,byrow = T),2,rev))))
  p8 <- ggplot2::ggplot(data=g, aes(x=Var1, y=Var2, fill=value))+ 
    geom_tile(color="white") +
    scale_fill_gradient2(low = "white", high = "red",midpoint=0.3,limit = c(0,1)) +
    ggtitle("HLA-C")
  p8
  
  
  
  top_gene_var <- gene_data[top_IDS_pod[1],]
  
  
  img <- do.call('cbind',lapply((c(0,1,2,3))^2,function(X){
    internal_top_gene_var <- top_gene_var
    internal_top_gene_var["SIX2"] <- internal_top_gene_var["SIX2"] + X
    internal_top_gene_var["NOTCH1"] <- internal_top_gene_var["NOTCH1"] + X
    internal_top_gene_var["WNT4"] <- internal_top_gene_var["WNT4"] + X
    
    internal_top_gene_var["FGF8"] <- internal_top_gene_var["FGF8"] + X
    internal_top_gene_var["FGFR1"] <- internal_top_gene_var["FGFR1"] + X
    
    c(t(internal_top_gene_var)%*%t(gene2pixel))
  }))
  
  g <- (convert_to_RGB(aperm(array((img[,1]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p9 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p9
  
  
  
  # img <- t(top_gene_var)%*%t(gene2pixel)
  g <- (convert_to_RGB(aperm(array((img[,2]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p10 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p10
  
  
  
  g <- (convert_to_RGB(aperm(array((img[,3]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p11 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p11
  
  
  g <- (convert_to_RGB(aperm(array((img[,4]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p12 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p12
  
  
  
  
  
  
  
  
  
  
  
  # top_gene_var <- gene_data[top_IDS_tub[64],]
  
  
  img <- do.call('cbind',lapply((c(0,3,6,9))^2,function(X){
    internal_top_gene_var <- top_gene_var
    internal_top_gene_var["CDH1"] <- internal_top_gene_var["CDH1"] + X
    internal_top_gene_var["PAX2"] <- internal_top_gene_var["PAX2"] + X
    internal_top_gene_var["GATA3"] <- internal_top_gene_var["GATA3"] + X
    
    
    c(t(internal_top_gene_var)%*%t(gene2pixel))
  }))
  
  g <- (convert_to_RGB(aperm(array((img[,1]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p13 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p13
  
  
  
  # img <- t(top_gene_var)%*%t(gene2pixel)
  g <- (convert_to_RGB(aperm(array((img[,2]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p14 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p14
  
  
  
  g <- (convert_to_RGB(aperm(array((img[,3]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p15 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p15
  
  
  g <- (convert_to_RGB(aperm(array((img[,4]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p16 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p16
  
  
  
  
  
  
  
  
  # top_gene_var <- gene_data[top_IDS_tub[64],]
  
  
  img <- do.call('cbind',lapply((c(0,3,6,9))^2,function(X){
    internal_top_gene_var <- top_gene_var
    
    internal_top_gene_var["IGKC"] <- internal_top_gene_var["IGKC"] + X
    internal_top_gene_var["CD74"] <- internal_top_gene_var["CD74"] + X
    internal_top_gene_var["B2M"] <- internal_top_gene_var["B2M"] + X
    
    c(t(internal_top_gene_var)%*%t(gene2pixel))
  }))
  
  g <- (convert_to_RGB(aperm(array((img[,1]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p17 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p17
  
  
  
  # img <- t(top_gene_var)%*%t(gene2pixel)
  g <- (convert_to_RGB(aperm(array((img[,2]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p18 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p18
  
  
  
  g <- (convert_to_RGB(aperm(array((img[,3]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p19 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p19
  
  
  g <- (convert_to_RGB(aperm(array((img[,4]),dim=c(3,32,32)),c(2,3,1))))
  g <- rasterGrob(g, interpolate=TRUE)
  p20 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
    annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
  p20
  
  
  
  
  
  gl <- list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20)
  
  gl <- lapply(gl,function(X){
    X + theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank())
    
  })
  
  lm <- rbind(c(17,13,9,1,1,1,2),
              c(18,14,10,1,1,1,3),
              c(19,15,11,1,1,1,4),
              c(20,16,12,5,6,7,8)
  )
  
  library(grid)
  library(gridExtra)
  gg1 <- arrangeGrob(
    grobs = gl,
    layout_matrix = lm
  )
  
  # ggsave(gg1,filename = "./figures/plot_showcase/showcase.png",width = 15,height = 8.5)
  
}


if (F){
  
  
  known_pixel_tiles <- pixel_nanostring.gcproc$transformed.data$x
  known_gene_tiles <- pixel_nanostring.gcproc$transformed.data$y

  known_pixel_tiles <- known_pixel_tiles[rowSums(pixel_nanostring)!=0,]
  known_gene_tiles <- known_gene_tiles[rowSums(pixel_nanostring)!=0,]
  
  
  v <- pixel_nanostring.gcproc$main.parameters$v.beta
  u <- pixel_nanostring.gcproc$main.parameters$u.beta
  
  set.seed(1)
  train_id <- sample(c(1:198),as.integer(198*0.8))
  test_id <- c(1:198)[-train_id]
  
  data.u.st_test <- known_pixel_tiles[test_id,]
  
  predicted_known_gene_tiles <- (data.u.st_test%*%t(pixel2gene*(pixel2gene>0)))
  predicted_known_gene_tiles <- predicted_known_gene_tiles + abs(min(predicted_known_gene_tiles))
  predicted_known_gene_tiles <- log2(1+predicted_known_gene_tiles*1e4 / rowSums(predicted_known_gene_tiles+1e-10))
  predicted_known_gene_tiles <- (predicted_known_gene_tiles*(mean(known_gene_tiles)/mean(predicted_known_gene_tiles)))
  colnames(predicted_known_gene_tiles) <- colnames(gene_data)
  
  cos.sim <- function(A,B) 
  {
    return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
  }   
  
  cor_cells <- do.call('c',lapply(c(1:length(test_id)),function(X){
    cos.sim((predicted_known_gene_tiles[X,]),(known_gene_tiles[test_id[X],]))
  }))
  
  cor_genes <- do.call('c',lapply(c(1:dim(predicted_known_gene_tiles)[2]),function(X){
    cos.sim((predicted_known_gene_tiles[,X]),(known_gene_tiles[test_id,X]))
  }))
  
  
  g1 <- ggplot(data.frame(Cosine_Correlation=cor_cells[!is.na(cor_cells)],All_Cells=" "), aes(x=Cosine_Correlation,y=All_Cells)) + geom_violin() 
  g2 <- ggplot(data.frame(Cosine_Correlation=cor_genes[!is.na(cor_genes)],All_Genes=" "), aes(x=Cosine_Correlation,y=All_Genes)) + geom_violin() 
  
}


if (F){
  
  all_list_ohe <- as.matrix(mltools::one_hot(data.table::data.table(data.frame(as.factor(all_list)))))
  
  v.sc <- unrelated_pixel_nanostring.gcproc$main.parameters$v.beta
  v <- pixel_nanostring.gcproc$main.parameters$v.beta
  u <- pixel_nanostring.gcproc$main.parameters$u.beta
  known_sc_gene_tiles <- unrelated_pixel_nanostring.gcproc$transformed.data$y
  
  
  beta_simil <- pixel_nanostring.gcproc$transformed.data$x%*%u%*%t(known_sc_gene_tiles%*%v)
  aligned_top_cell_types <- do.call('c',lapply(c(1:dim(beta_simil)[1]),function(X){main_kidney_list[all_list[which(beta_simil[X,]==max(beta_simil[X,]))[1]],2]}))[grep("PanCK",sort(row.names(pixel_nanostring.gcproc[["transformed.data"]][["y"]])))]

  g0<-ggplot2::ggplot(data=data.frame(Observed_Renal_Tubule_Spots = c(table(aligned_top_cell_types)), Predicted_Cell_Types = names(table(aligned_top_cell_types))), aes(x=Observed_Renal_Tubule_Spots, y=Predicted_Cell_Types)) +
    geom_bar(stat="identity")
  
  
  glm <- rbind(c(1),
              c(1),
              c(2),
              c(3)
  )
  
  library(grid)
  library(gridExtra)
  gg2 <- arrangeGrob(
    grobs = list(g0,g1,g2),
    layout_matrix = lm
  )
  
  # ggsave(gg2,filename = "./figures/plot_metrics/metrics.png",width = 5,height = 5)
  
  
}




if (F){
  
  
  convert_to_COMP <- function(x){
    e <- ecdf(x)
    j <- e(x)
    x <- array(j,dim=dim(x))
    return(x)
  }
  
  library(ggplot2)
  library(grid)
  plot_list <- lapply(c(1:6),function(i){
    component <- unrelated_pixel_nanostring.gcproc$main.parameters$u.beta[,i]
    
    g <- (convert_to_COMP(aperm(array((component),dim=c(3,32,32)),c(2,3,1))))
    g <- rasterGrob(g, interpolate=TRUE)
    p9 <- ggplot2::qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) 
    p9
  })
  
  
  plot_list <- lapply(plot_list,function(X){
    X + theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank())
    
  })
  
  
  lm <- rbind(c(1,4),
              c(2,5),
              c(3,6)
  )
  
  library(grid)
  library(gridExtra)
  gg <- arrangeGrob(
    grobs = plot_list,
    layout_matrix = lm
  )
  
  # ggsave(gg,filename = "./figures/plot_workflow/main_components.png")
  
}

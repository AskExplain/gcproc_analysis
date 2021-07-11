gene_data <- read.csv("./data/workflow/main_image_as_gene_v1.csv",header=T)
cell_labels <- read.csv("./data/workflow/main_kidney_list.csv",header=T)
gene_type <- gene_data[,1]
cell_labels <- cell_labels[gene_type,2]
gene_data <- t(gene_data[,-1])
gene2pixel <- read.csv("./data/workflow/gene2pixel.csv",header=T,row.names=1)


domain <- c(1:50)


set.seed(1)

main_all <- do.call('cbind',lapply(c(1:5),function(X){
  L <- runif(1,30,50)
  x0 <- runif(1,10,40)
  k <- runif(1,0.15,0.35)
  val <- L/(1+exp(-k*(domain-x0)))
  plot(val)
  return(val)
}))





gene_id=1
main_gene <- gene_data[,gene_id]
names(main_gene) <- row.names(gene_data)
genes <- c("SIX2","NPHS2","WT1","PODXL","NPHS1")


convert_to_RGB <- function(X){
  e <- ecdf(X)
  j <- e(X)
  X <- array(j,dim=dim(X))
  return(X)
}

library(gridExtra)
library(grid)
library(ggplot2)

for (i in 1:max(domain)){
  print(i)
  main_gene <- gene_data[,gene_id]
  
  print(i)
  internal_gene <- main_gene
  internal_gene[genes[1]] <- main_gene[genes[1]] + main_all[i,1]
  internal_gene[genes[2]] <- main_gene[genes[2]] + main_all[i,2]
  internal_gene[genes[3]] <- main_gene[genes[3]] + main_all[i,3]
  internal_gene[genes[4]] <- main_gene[genes[4]] + main_all[i,4]
  internal_gene[genes[5]] <- main_gene[genes[5]] + main_all[i,5]
  
  cell_image <- t(internal_gene)%*%t(gene2pixel)
  
  final_cell_image <- aperm(array(convert_to_RGB(cell_image),dim=c(3,32,32)),c(2,3,1))
  
  
  g <- rasterGrob(final_cell_image, interpolate=TRUE)
  
  gs <- lapply(c(1:dim(main_all)[2],6),function(X){
    if (X==6){
      a <- qplot(1:10, 1:10, geom="blank") +
      annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
    } else {
      a <- ggplot(data = data.frame(x=domain[1:i],y=main_gene[genes[X]] + main_all[1:i,X]),
                  aes(x=x, y=y)) +
        geom_line()+
        xlim(0, max(domain))+
        ylim(0, 100)+
        ggtitle(genes[X])
    }
    return(a)
  })
  
  lay <- rbind(c(1,6,6,6),
               c(2,6,6,6),
               c(3,6,6,6),
               c(4,6,6,6),
               c(5,6,6,6))
  
  p <- arrangeGrob(grobs = gs, layout_matrix = lay)

  ggsave(file=paste("./figures/visual_gene_regulatory_network/main___",i,".png",sep=""),p,width=14,height=10)
  
  print(cell_labels[gene_id])
  }



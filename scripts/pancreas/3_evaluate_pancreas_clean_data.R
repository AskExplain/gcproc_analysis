# Loading in the functions.

# Sourced from MNN paper and github
# https://github.com/MarioniLab/FurtherMNN2018/blob/master/simulations/functions.R
source("functions.R")

save("./data/workflow/run_pancreas_clean_data.Rdata")

out1 <- runAllMethods(b, x)
out2 <- runAllMethods(b, s)

out3 <- lapply(c(1:5),function(X){
  data <- rbind(out1$mat[[X]],out2$mat[[X]][-c(1:dim(b)[2]),])
  umap_data <- umap::umap(data)$layout
  row.names(umap_data) <- row.names(data)
  return(list(data=data,umap_data=umap_data))
})

names(out3) <- names(out1$mat)
out3$gencoder$data <- rbind(b.b%*%(v*u*w),x.x%*%(v*u*w),s.s%*%(v*u*w))
out3$gencoder$umap_data <- umap::umap(rbind(b.b%*%(v*u*w),x.x%*%(v*u*w),s.s%*%(v*u*w)))$layout

names(out3)[6] <- "Generative Encoder"


library(ggplot2)

gl <- lapply(c(1:6),function(X){
  gg <- ggplot2::ggplot(data.frame(V1=out3[[X]][["umap_data"]][,1],V2=out3[[X]][["umap_data"]][,2],cell_types = row.names(out3[[6]][["umap_data"]])), aes(x=V1, y=V2,col=cell_types)) +
    geom_point(size=0.1, shape=19) +
    ggtitle(names(out3)[X])
  
  if (X!=3){
    gg <- gg + 
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),legend.position="none",
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank())
    
  } else {
    
    gg <- gg  + 
      theme(axis.line=element_blank(),axis.text.x=element_blank(),
            axis.text.y=element_blank(),axis.ticks=element_blank(),
            axis.title.x=element_blank(),
            axis.title.y=element_blank(),
            panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),plot.background=element_blank())
    
  }
  return(gg)
})


lm <- rbind(c(1,2,3),c(4,5,6))

library(grid)
library(gridExtra)
gg1 <- arrangeGrob(
  grobs = gl,
  layout_matrix = lm
)

ggsave(gg1,filename = "../../../../figures/plot_benchmark/pancreas.png",width = 8,height = 4)


















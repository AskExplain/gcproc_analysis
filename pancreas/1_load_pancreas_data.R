num_genes = 2000
seed = 1

set.seed(seed)

setwd("")

baron <- readRDS("./data/workflow/main_rds/baron.rds")
baron_gene_order <- order(apply(baron@assays$data$counts,1,var),decreasing=T)
baron_gene_IDS <- baron_gene_order[1:num_genes]
baron_gene_names <- row.names(baron)
baron_gene_names_IDS <- baron_gene_names[baron_gene_IDS]


segerstolpe <- readRDS("./data/workflow/main_rds/segerstolpe.rds")
segerstolpe_gene_order <- order(apply(segerstolpe@assays$data$counts,1,var),decreasing=T)
segerstolpe_gene_IDS <- segerstolpe_gene_order[1:num_genes]
segerstolpe_gene_names <- row.names(segerstolpe)


xin <- readRDS("./data/workflow/main_rds/xin.rds")
xin_gene_order <- order(apply(xin@assays$data$normcounts,1,var),decreasing=T)
xin_gene_IDS <- xin_gene_order[1:num_genes]
xin_gene_names <- row.names(xin)


top_gene_intersect <- Reduce('intersect',list(baron_gene_names[baron_gene_order],segerstolpe_gene_names[segerstolpe_gene_order],xin_gene_names[xin_gene_order]))[1:num_genes]

baron_gene_IDS <- which(baron_gene_names %in% top_gene_intersect)
segerstolpe_gene_IDS <- which(segerstolpe_gene_names %in% top_gene_intersect)
xin_gene_IDS <- which(xin_gene_names %in% top_gene_intersect)

# Baron
clusters <- as.factor(baron$cell_type1)
cluster_list <- c("alpha","beta","delta","gamma","acinar","epsilon","endothelial")
counts <- baron@assays$data$counts[,clusters%in%c(cluster_list)]
baron_clusters <- as.character(clusters[clusters%in%c(cluster_list)])

holds_covariates = T

cluster_order <- order(unique(baron_clusters))

baron_data <-
  lapply(c(1:length(unique(baron_clusters))),function(Y){
    
    ref_y <- counts[baron_gene_IDS,baron_clusters == unique(baron_clusters)[cluster_order[Y]]]
    order_ref_y <- order(apply(ref_y,2,var),decreasing = T)
    colnames(ref_y) <- baron_clusters[baron_clusters == unique(baron_clusters)[cluster_order[Y]]]
    
    
    ref_y <- (ref_y[,order_ref_y])
    
    row.names(ref_y) <- c(sort(baron_gene_names[baron_gene_IDS]))
    
    return(ref_y)
  })

names(baron_data) <- sort(unique(baron_clusters))



Z_baron_data <-
  lapply(c(1:length(unique(baron_clusters))),function(Y){
    
    ref_y <- counts[baron_gene_IDS,baron_clusters == unique(baron_clusters)[cluster_order[Y]]]
    order_ref_y <- order(apply(ref_y,2,var),decreasing = T)
    colnames(ref_y) <- baron_clusters[baron_clusters == unique(baron_clusters)[cluster_order[Y]]]
    
    
    ref_y <- (ref_y[,order_ref_y])
    
    row.names(ref_y) <- c(sort(baron_gene_names[baron_gene_IDS]))
    
    return(ref_y)
  })

names(Z_baron_data) <- sort(unique(baron_clusters))



# Segerstolpe
clusters <- as.factor(segerstolpe$cell_type1)
new_clusters <- segerstolpe@colData$disease

cluster_list <- c("alpha","beta","delta","gamma","acinar","epsilon","endothelial")
new_cluster_list <- c("normal")

covariates <- as.data.frame(segerstolpe@colData[,c(4,5)])
covariates$sex <- as.factor(covariates$sex)
covariates <- mltools::one_hot(data.table::as.data.table(covariates))

covariates <- covariates[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list),]
s.c <- covariates
counts <- segerstolpe@assays$data$counts[,clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)]
segerstolpe_clusters <- as.character(clusters[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)])

holds_covariates = T

cluster_order <- order(unique(segerstolpe_clusters))

Z_segerstolpe_data <-
  lapply(c(1:length(unique(segerstolpe_clusters))),function(Y){
    
    ref_y <- counts[sort(segerstolpe_gene_names[segerstolpe_gene_IDS]),segerstolpe_clusters == unique(segerstolpe_clusters)[cluster_order[Y]]]
    order_ref_y <- order(apply(ref_y,2,var),decreasing = T)
    colnames(ref_y) <- segerstolpe_clusters[segerstolpe_clusters == unique(segerstolpe_clusters)[cluster_order[Y]]]
    
    if (holds_covariates == T){
      ref_y_covariates <- t(covariates[segerstolpe_clusters == unique(segerstolpe_clusters)[cluster_order[Y]],])
      ref_y_covariates <- as.matrix(rbind(ref_y_covariates)[,order_ref_y])
      ref_y <- as.matrix(rbind(ref_y)[,order_ref_y])
      
    }
    else {
      
      ref_y <- (ref_y[,order_ref_y])
      
      row.names(ref_y) <- c(sort(segerstolpe_gene_names[segerstolpe_gene_IDS]))
      
    }
    
    return(list(ref_y=ref_y,ref_y_covariates=ref_y_covariates))
  })

names(Z_segerstolpe_data) <- sort(unique(segerstolpe_clusters))



# Xin
clusters <- as.factor(xin$cell_type1)
new_clusters <- as.factor(xin@colData$condition)

cluster_list <- c("alpha","beta","delta","gamma","acinar","epsilon","endothelial")
new_cluster_list <- c("Healthy")

covariates <- as.data.frame(xin@colData[,-c(1,2,6)])
covariates$gender <- as.factor(covariates$gender)
covariates$ethnicity <- as.factor(covariates$ethnicity)
covariates <- mltools::one_hot(data.table::as.data.table(covariates))

covariates <- covariates[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list),]
x.c <- covariates
counts <- xin@assays$data$normcounts[,clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)]
xin_clusters <- as.character(clusters[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)])

holds_covariates = T

cluster_order <- order(unique(xin_clusters))

Z_xin_data <-
  lapply(c(1:length(unique(xin_clusters))),function(Y){
    
    ref_y <- counts[sort(xin_gene_names[xin_gene_IDS]),xin_clusters == unique(xin_clusters)[cluster_order[Y]]]
    order_ref_y <- order(apply(ref_y,2,var),decreasing = T)
    colnames(ref_y) <- xin_clusters[xin_clusters == unique(xin_clusters)[cluster_order[Y]]]
    
    if (holds_covariates == T){
      ref_y_covariates <- t(covariates[xin_clusters == unique(xin_clusters)[cluster_order[Y]],])
      ref_y_covariates <- as.matrix(rbind(ref_y_covariates)[,order_ref_y])
      ref_y <- as.matrix(rbind(ref_y)[,order_ref_y])
      
    }
    else {
      
      ref_y <- (ref_y[,order_ref_y])
      
      row.names(ref_y) <- c(sort(xin_gene_names[xin_gene_IDS]))
      
    }
    
    return(list(ref_y=ref_y,ref_y_covariates=ref_y_covariates))
  })

names(Z_xin_data) <- sort(unique(xin_clusters))

b <- (do.call('cbind',lapply(baron_data,function(X){X})))
s <- (do.call('cbind',lapply(Z_segerstolpe_data,function(X){X$ref_y})))
x <- (do.call('cbind',lapply(Z_xin_data,function(X){X$ref_y})))

b <- log2(1+b*10000/colSums(b))
x <- log2(1+x*10000/colSums(x))
s <- log2(1+s*10000/colSums(s))


save(list = c("b","x","x.c","s","s.c"),file="./data/workflow/load_pancreas_clean_data_all.Rdata")





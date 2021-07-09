num_genes = 2000
seed = 1

set.seed(seed)


baron <- readRDS("./data/workspace/main_rds/baron.rds")
baron_gene_order <- order(apply(baron@assays$data$counts,1,var),decreasing=T)
baron_gene_IDS <- baron_gene_order[1:num_genes]
baron_gene_names <- row.names(baron)
baron_gene_names_IDS <- baron_gene_names[baron_gene_IDS]


segerstolpe <- readRDS("./data/workspace/main_rds/segerstolpe.rds")
segerstolpe_gene_order <- order(apply(segerstolpe@assays$data$counts,1,var),decreasing=T)
segerstolpe_gene_IDS <- segerstolpe_gene_order[1:num_genes]
segerstolpe_gene_names <- row.names(segerstolpe)


xin <- readRDS("./data/workspace/main_rds/xin.rds")
xin_gene_order <- order(apply(xin@assays$data$normcounts,1,var),decreasing=T)
xin_gene_IDS <- xin_gene_order[1:num_genes]
xin_gene_names <- row.names(xin)


top_gene_intersect <- Reduce('intersect',list(baron_gene_names[baron_gene_order],segerstolpe_gene_names[segerstolpe_gene_order],xin_gene_names[xin_gene_order]))[1:num_genes]

baron_gene_IDS <- which(baron_gene_names %in% top_gene_intersect)
segerstolpe_gene_IDS <- which(segerstolpe_gene_names %in% top_gene_intersect)
xin_gene_IDS <- which(xin_gene_names %in% top_gene_intersect)

# Baron

# genes_interest <- c("KIAA1324","GCGR","INSR")
# gene_IDS <- unique(c(gene_IDS,which(gene_names%in%genes_interest)))


clusters <- as.factor(baron$cell_type1)
cluster_list <- c("alpha","beta","delta","gamma")
counts <- baron@assays$data$counts[,clusters%in%c(cluster_list)]
baron_clusters <- as.character(clusters[clusters%in%c(cluster_list)])
baron_data <- counts[baron_gene_IDS,]


# Segerstolpe
clusters <- as.factor(segerstolpe$cell_type1)
new_clusters <- segerstolpe@colData$disease

cluster_list <- c("alpha","beta","delta","gamma")
new_cluster_list <- c("normal")

covariates <- as.data.frame(segerstolpe@colData[,c(4,5)])
covariates$sex <- as.factor(covariates$sex)
covariates <- mltools::one_hot(data.table::as.data.table(covariates))

segerstolpe_covariates <- covariates[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list),]
counts <- segerstolpe@assays$data$counts[,clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)]
segerstolpe_clusters <- as.character(clusters[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)])

segerstolpe_data <- counts[segerstolpe_gene_IDS,]



# Xin

# genes_interest <- c("KIAA1324","GCGR","INSR")
# gene_IDS <- unique(c(gene_IDS,which(gene_names%in%genes_interest)))

clusters <- as.factor(xin$cell_type1)
new_clusters <- as.factor(xin@colData$condition)

cluster_list <- c("alpha","beta","delta","gamma")
new_cluster_list <- c("Healthy")

covariates <- as.data.frame(xin@colData[,-c(1,2,6)])
covariates$gender <- as.factor(covariates$gender)
covariates$ethnicity <- as.factor(covariates$ethnicity)
covariates <- mltools::one_hot(data.table::as.data.table(covariates))

xin_covariates <- covariates[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list),]
counts <- xin@assays$data$normcounts[,clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)]
xin_clusters <- as.character(clusters[clusters%in%c(cluster_list) & new_clusters%in%c(new_cluster_list)])
xin_data <- counts[xin_gene_IDS,]

b <- baron_data
s <- segerstolpe_covariates
s.c <- segerstolpe_data

x <- xin_data
x.c <- xin_covariates

b <- log2(1+b*10000/colSums(b))
x <- log2(1+x*10000/colSums(x))
s <- log2(1+s*10000/colSums(s))

x.c_fix <- x-(x%*%t(x.c)%*%MASS::ginv(x.c%*%t(x.c))%*%x.c)
s.c_fis <- s-(s%*%t(s.c)%*%MASS::ginv(s.c%*%t(s.c))%*%s.c)

save(list = c("b","x.c_fix","s.c_fix"),"./data/workflow/load_pancreas_clean_data.Rdata")





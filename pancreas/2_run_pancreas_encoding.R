load("./data/workflow/load_pancreas_clean_data_all.Rdata")

# Load gcproc
# library(gcproc)
# Alternatively, load via devtools
# devtools::load_all("./gcproc/R/")

print("Frequency of cell types per dataset")
print(table(data.frame(col=as.factor(colnames(cbind(b,x,s))),pch=c(rep(1,dim(b)[2]),rep(2,dim(x)[2]),rep(3,dim(s)[2])))))


gcproc.config <- list(k_dim = 70,
                      j_dim = 70,
                      eta=1e-2,
                      max_iter=2500,
                      min_iter = 25,
                      tol=1,
                      log=F,
                      center=T,
                      scale.z=T,
                      batches=16,
                      cores=8,
                      verbose=T,
                      init="svd-quick")


x.s_cv.gcproc_cov <- gcproc::gcproc(
  y = t(x),
  x = t(s),
  config = gcproc.config,
  seed = seed,
  covariates = list(covariates_y.sample = t(x.c),
                    covariates_y.feature = NULL,
                    covariates_x.sample = t(s.c),
                    covariates_x.feature = NULL
  )
  
)


x.b_cv.gcproc_cov <- gcproc::gcproc(
  y = t(s),
  x = t(b),
  seed = seed,
  config = gcproc.config,
  covariates = list(covariates_y.sample = t(s.c),
                    covariates_y.feature = NULL,
                    covariates_x.sample = NULL,
                    covariates_x.feature = NULL
  )
)



b.b <- t(b)
x.x <- t(x) 
s.s <- t(s) 

v <- x.s_cv.gcproc_cov$main.parameters$v.beta
u <- x.s_cv.gcproc_cov$main.parameters$u.beta
w <- x.b_cv.gcproc_cov$main.parameters$u.beta
r <- x.b_cv.gcproc_cov$main.parameters$v.beta

plot(umap::umap(x.x%*%(v*u+w*r))$layout,col=as.factor(row.names(x.x)))
plot(umap::umap(s.s%*%(v*u+w*r))$layout,col=as.factor(row.names(s.s)))
plot(umap::umap(b.b%*%(v*u+w*r))$layout,col=as.factor(row.names(b.b)))

plot(umap::umap(rbind(b.b,x.x,s.s)%*%(v*u+w*r))$layout,col=as.factor(row.names(rbind(b.b,x.x,s.s))),pch=c(rep(1,dim(b.b)[1]),rep(2,dim(x.x)[1]),rep(3,dim(s.s)[1])),cex=0.5)

save(list=c("u","v","w","r","b.b","x.x","s.s"),file="./data/workflow/run_pancreas_clean_data.Rdata")

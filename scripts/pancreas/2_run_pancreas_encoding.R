load("./data/workflow/load_pancreas_clean_data.Rdata")

library(gcproc)

print("Frequency of cell types per dataset")
print(table(data.frame(col=as.factor(colnames(cbind(b,x,s))),pch=c(rep(1,dim(b)[2]),rep(2,dim(x)[2]),rep(3,dim(s)[2])))))


gcproc.config <- list(k_dim = 70,
                      j_dim = 70,
                      eta=5e-3,
                      max_iter=1500,
                      min_iter = 15,
                      tol=1e-8,
                      log=F,
                      center=T,
                      scale.z=T,
                      batches=128,
                      cores=8,
                      verbose=T,
                      init="svd")




x.s_cv.gcproc_cov <- gcproc::cv.gcproc(
  y = t(x),
  x = t(s),
  config = gcproc.config,
  initial_starts = 10
    )

x.b_cv.gcproc_cov <- gcproc::cv.gcproc(
  y = t(x),
  x = t(b),
  initial_starts = 10,
  config = gcproc.config,
  anchors = list(  anchor_y.sample = NULL,
                   anchor_y.feature = x.b_cv.gcproc_cov$main.parameters$v.beta,
                   anchor_x.sample = NULL,
                   anchor_x.feature = NULL
                )
)




b.b <- x.b_cv.gcproc_cov$transformed.data$x
x.x <- x.s_cv.gcproc_cov$transformed.data$y
s.s <- x.s_cv.gcproc_cov$transformed.data$x

v <- x.s_cv.gcproc_cov$main.parameters$v.beta
u <- x.s_cv.gcproc_cov$main.parameters$u.beta
w <- x.b_cv.gcproc_cov$main.parameters$u.beta

plot(umap::umap(b.b%*%(u*v*w))$layout,col=as.factor(row.names(b.b)))
plot(umap::umap(x.x%*%(u*w*v))$layout,col=as.factor(row.names(x.x)))
plot(umap::umap(s.s%*%(u*w*v))$layout,col=as.factor(row.names(s.s)))

plot(umap::umap(rbind(b.b,x.x,s.s)%*%(u*v*w))$layout,col=as.factor(row.names(rbind(b.b,x.x,s.s))),pch=c(rep(1,dim(b.b)[1]),rep(2,dim(x.x)[1]),rep(3,dim(s.s)[1])),cex=0.5)

save(list=c("u","v","w","b.b","x.x","s.s"),file="./data/workflow/run_pancreas_clean_data.Rdata")

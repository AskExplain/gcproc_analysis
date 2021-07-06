load("./data/workflow/load_pancreas_clean_data.Rdata")
devtools::load_all("~/Documents/main_files/GCP/explain-dpm/clean/gcproc_generalised_canonical_procrustes_2020_06_04/gcproc//R/")

print("Frequency of cell types per dataset")
print(table(data.frame(col=as.factor(row.names(rbind(b,x.c_fix,s.c_fix))),pch=c(rep(1,dim(b.b)[1]),rep(2,dim(x.c_fix)[1]),rep(3,dim(s.c_fix)[1])))))


gcproc.config <- list(k_dim = 70,
                      j_dim = 70,
                      eta=1e-1,
                      max_iter=1500,
                      min_iter = 15,
                      tol=1e-8,
                      log=F,
                      center=T,
                      scale.z=T,
                      batches=64,
                      cores=8,
                      verbose=T,
                      init="svd")


x.x_cv.gcproc_cov <- gcproc::cv.gcproc(
  y = t(x-(x%*%t(x.c)%*%MASS::ginv(x.c%*%t(x.c))%*%x.c)),
  x = t(x-(x%*%t(x.c)%*%MASS::ginv(x.c%*%t(x.c))%*%x.c)),
  initial_starts = 1,
  config = gcproc.config
)


s.s_cv.gcproc_cov <- gcproc::cv.gcproc(
  y = t(s-(s%*%t(s.c)%*%MASS::ginv(s.c%*%t(s.c))%*%s.c)),
  x = t(s-(s%*%t(s.c)%*%MASS::ginv(s.c%*%t(s.c))%*%s.c)),
  initial_starts = 1,
  config = gcproc.config
)

x.s.config <- gcproc.config
x.s.config$tol = 1e-3

x.s_cv.gcproc_cov <- gcproc::transfer.gcproc(
  y = t(x-(x%*%t(x.c)%*%MASS::ginv(x.c%*%t(x.c))%*%x.c)),
  x = t(s-(s%*%t(s.c)%*%MASS::ginv(s.c%*%t(s.c))%*%s.c)),
  initial_starts = 1,
  config = x.s.config,
  anchors = list( anchor_y.sample = x.x_cv.gcproc_cov$main.parameters$alpha.L.J*x.x_cv.gcproc_cov$main.parameters$alpha.L.K,
                  anchor_y.feature = NULL,
                  anchor_x.sample = s.s_cv.gcproc_cov$main.parameters$alpha.L.J*s.s_cv.gcproc_cov$main.parameters$alpha.L.K,
                  anchor_x.feature = NULL  )
)


b.config <- list(k_dim = 70,
                 j_dim = 70,
                 eta=1e-1,
                 max_iter=0,
                 min_iter = 0,
                 tol=1e-3,
                 log=F,
                 center=T,
                 scale.z=T,
                 batches=64,
                 cores=8,
                 verbose=T,
                 init="random")


b.b_cv.gcproc_cov <- gcproc::cv.gcproc(
  y = t(b),
  x = t(b),
  initial_starts = 1,
  config = gcproc.config
)

x.b_cv.gcproc_cov <- gcproc::transfer.gcproc(
  y = t(x-(x%*%t(x.c)%*%MASS::ginv(x.c%*%t(x.c))%*%x.c)),
  x = t(b),
  config = gcproc.config,
  initial_starts = 1,
  anchors = list( anchor_y.sample = x.x_cv.gcproc_cov$main.parameters$alpha.L.J*x.x_cv.gcproc_cov$main.parameters$alpha.L.K,
                  anchor_y.feature = x.s_cv.gcproc_cov$main.parameters$v.beta*x.s_cv.gcproc_cov$main.parameters$u.beta,
                  anchor_x.sample = b.b_cv.gcproc_cov$main.parameters$alpha.L.J*b.b_cv.gcproc_cov$main.parameters$alpha.L.K,
                  anchor_x.feature = NULL  )
)



model.gcp1 <- x.s_cv.gcproc_cov

b.b <- b.b_cv.gcproc_cov$transformed.data$x
x.x <- model.gcp1$transformed.data$y
s.s <- model.gcp1$transformed.data$x

v <- model.gcp1$main.parameters$v.beta
u <- model.gcp1$main.parameters$u.beta

model.gcp2 <- x.b_cv.gcproc_cov
w <- model.gcp2$final_rotation

plot(umap::umap(b.b%*%(u*v*w))$layout,col=as.factor(row.names(b.b)))
plot(umap::umap(x.x%*%(u*v*w))$layout,col=as.factor(row.names(x.x)))
plot(umap::umap(s.s%*%(u*v*w))$layout,col=as.factor(row.names(s.s)))

plot(umap::umap(rbind(rbind(b.b),x.x,s.s)%*%(u*v*w))$layout,col=as.factor(row.names(rbind(b.b,x.x,s.s))),pch=c(rep(1,dim(b.b)[1]),rep(2,dim(x.x)[1]),rep(3,dim(s.s)[1])),cex=0.5)

x <- x.c_fix
s <- s.c_fix

save(list=c("u","v","w","b.b","x.x","s.s","b","x"),file="./data/workflow/run_pancreas_clean_data.Rdata")

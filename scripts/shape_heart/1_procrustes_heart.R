# Simulation sourced from
# https://rstudio-pubs-static.s3.amazonaws.com/222019_a2c617e90530422285420091c090d415.html


# number of data points
k = 500

# x
set.seed(7)
x <- rnorm(n = k, mean = 0, sd = 1)

# random noise
set.seed(9)
e <- rnorm(n = k, mean = 0, sd = 1)

# y
y <- abs(x) + e

# Matrix A
# Store x and y into a 2 x k matrix
A <- matrix(c(x, y), byrow = TRUE, nrow = 2)

# matrix B
# Arbitrarily rotate, scale, and translate A 
theta <- pi / 4
rot.mtx <- matrix(c(sin(theta), -cos(theta), cos(theta), sin(theta)), ncol=2)
B <- 0.5*rot.mtx %*% A - 6


devtools::load_all("~/Documents/main_files/GCP/explain-dpm/clean/gcproc_generalised_canonical_procrustes_2020_06_04/gcproc/R/")
A.A_cv.gcproc <- gcproc::gcproc(
  y = t(A),
  x = t(A),
  seed = 1,
  config = list(k_dim = 15, j_dim = 2, eta = 1e-2, max_iter = 1500, min_iter = 15,
                tol = 0.01, log = F, center = T, scale.z = T, batches = 1, cores = 1, verbose = T,
                init = "eigen-dense")
)

B.B_cv.gcproc <- gcproc::gcproc(
  y = t(B),
  x = t(B),
  seed = 1,
  config = list(k_dim = 15, j_dim = 2, eta = 1e-2, max_iter = 1500, min_iter = 15,
                tol = 0.01, log = F, center = T, scale.z = T, batches = 1, cores = 1, verbose = T,
                init = "eigen-dense")
)


A.B_cv.gcproc <- gcproc::gcproc(
  y = t(A),
  x = t(B),
  seed = 1,
  config = list(k_dim = 15, j_dim = 2, eta = 1e-2, max_iter = 1500, min_iter = 15,
                tol = 0.01, log = F, center = T, scale.z = T, batches = 1, cores = 1, verbose = T,
                init = "eigen-dense"),
  anchors = list(
    anchor_y.sample = A.A_cv.gcproc$main.parameters$alpha.L.J*A.A_cv.gcproc$main.parameters$alpha.L.K,
    anchor_y.feature = NULL,
    anchor_x.sample = B.B_cv.gcproc$main.parameters$alpha.L.J*B.B_cv.gcproc$main.parameters$alpha.L.K,
    anchor_x.feature = B.B_cv.gcproc$main.parameters$u.beta*B.B_cv.gcproc$main.parameters$v.beta)
)




procrustes <- function(A, B){
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}

procrustes.results <- procrustes(A, B)

test.Y <- A.B_cv.gcproc$transformed.data$y%*%(A.B_cv.gcproc$main.parameters$v.beta)
test.R.X <- A.B_cv.gcproc$transformed.data$x%*%A.B_cv.gcproc$main.parameters$u.beta + (A.B_cv.gcproc$transformed.data$y%*%A.B_cv.gcproc$main.parameters$v.beta - A.B_cv.gcproc$transformed.data$x%*%A.B_cv.gcproc$main.parameters$u.beta)

A.df <- as.data.frame(t(A))
B.df <- as.data.frame(t(B))

A.normalized.df <- (test.Y)
B.transformed.df <- (test.R.X)

A.pc.df <- (t(procrustes.results$A.normalized))
B.pc.df <- (t(procrustes.results$B.transformed))


data.df <- rbind(A.df, B.df, A.normalized.df, B.transformed.df, A.pc.df, B.pc.df)
colnames(data.df) <- c('x', 'y')
data.df$matrix <- rep(c('A', 'B', 'A', 'B', 'A', 'B'), each = k)
data.df$treatment <- rep(c('1_Original', '2_Generative Encoding', '3_Procrustes'), each = 2*k)

library(ggplot2)
data.plot <- ggplot(data.df, aes(x = x, y = y, color = matrix)) +
  geom_point(aes(shape = matrix), size = 4, alpha = 0.5) +
  theme_bw() +
  theme(aspect.ratio = 1) +
  facet_wrap(~treatment, scales = 'free') +
  labs(title = 'Procrustes analysis')

data.plot
# ggsave("~/Documents/main_files/Explain/gencoder/project_gencoder/figures/Figure_2_Statistical_Model_for_Generative_Encoding.png",data.plot,width = 7,height = 3)



  
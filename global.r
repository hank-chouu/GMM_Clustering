
library(MASS)
library(mclust)




# example data

precip <- data.frame(precip=log(read.csv('precip.csv')$PrecpMonth))
two_dim_example <- read.csv('2_dim_example.csv')
three_dim_example <- read.csv('three_dim_example.csv')

# gmm em functions

para_initialize <- function(X, k, method){ # X should be a matrix
  
  d <- dim(X)[2]
  w <- rep(1/k, k)
  
  init_mu <- matrix(0, d, k)
  
  if (method == 'runif'){
    
    for (i in 1:k){
      for (j in 1:d){
        init_mu[j, i] = runif(1, min(X[, j]), max(X[, j]))
      }
    }
  }
  
  if (method == 'kmeans'){
    init_mu <- t(kmeans(X, centers = k)$centers)
  }
  
  init_sigma <- list()
  for (i in 1:k){
    init_sigma[[i]] <- diag(d)
  }
  para <- list(weight = w, mu = init_mu, sigma = init_sigma)
  return(para)
}

# the form of para
# weight: k components; mu: row is dim, col is components; 
# sigma: list of k covariance matrices


e_step <- function(X, para, k){
  
  n <- dim(X)[1]
  z <-  matrix(NA, nrow=n, ncol= k)
  sum_of_comps <- rep(0, n)
  
  for (i in 1:k){
    dmv <- dmvnorm(X, para$mu[, i], para$sigma[[i]])
    comp <- dmv * para$weight[i]
    sum_of_comps <- sum_of_comps + comp
  }
  
  for (i in 1:k){
    dmv <- dmvnorm(X, para$mu[, i], para$sigma[[i]])
    comp_post <- dmv * para$weight[i] / sum_of_comps
    z[, i] <- comp_post
  }
  
  sum_of_comps_ln <- log(sum_of_comps)
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)
  
  list(z=z, lld=sum_of_comps_ln_sum)
}

m_step <- function(X, z, k){
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  w.next <- rep(0, k)
  mu.next <- matrix(0, d, k)
  sigma.next <- list()
  
  for (i in 1:k){
    
    new.w <- sum(z[, i])/n
    w.next[i] <- new.w
    
    new.mu <- 1/sum(z[, i])*colSums(z[, i]* X)
    mu.next[, i] <- new.mu
    
    new.mu <- matrix(rep(new.mu, n), n, d, byrow = T)
    
    new.sigma <- 1/sum(z[, i]) * t(z[, i]*(X - new.mu)) %*% (z[, i]*(X - new.mu))
    sigma.next[[i]] <- new.sigma
  }
  
  para.next <- list(weight = w.next, mu = mu.next, sigma = sigma.next)
  return(para.next)
}


generalize_gmm <- function(x, k, init.method){
  
  X <- as.matrix(x)
  para <- para_initialize(X, k, init.method)
  
  e.step <- e_step(X, para, k) # return z, lld
  para <- m_step(X, e.step$z, k) # return new para
  cur.loglik <- e.step$lld
  loglik.vector <- e.step$lld
  loglik.diff <- 1
  count <- 1
  
  while(loglik.diff > 1e-6){
    
    e.step <- e_step(X, para, k)
    para <- m_step(X, e.step$z, k)
    
    count <- count + 1
    # print(count)
    if (count > 300){
      stop('reached max iters')
    }
    loglik.vector <- c(loglik.vector, cur.loglik)
    loglik.diff <- abs((cur.loglik - e.step$lld))
    cur.loglik <- e.step$lld
    
  }
  cluster <- factor(apply(e.step$z, 1, which.max))
  list(para=para, z=e.step$z, count=count, lld_vector=loglik.vector, lld=cur.loglik, cluster=cluster, success=T)
}


wrapped_gmm <- function(X, k){
  
  count <- 0
  max_count <- 10
  
  while (count <= max_count) {
    res <- try({
      result <- generalize_gmm(X, k, 'kmeans')
      result
    }, silent = TRUE)
    count <- count + 1
    if (!inherits(res, "try-error")) {
      return(res)
      break
    } 
  }
  print('go for runif')
  count <- 0
  while (count <= 300) {
    res <- try({
      result <- generalize_gmm(X, k, 'runif')
      result
    }, silent = TRUE)
    count <- count + 1
    if (!inherits(res, "try-error")) {
      return(res)
      break
    } 
  }
  # stop('unable to cluster')
  list(para='unable to cluster', success=F)
}

#####--------------------------------------------
##### preamble

rm(list=ls())

library(gaussquad)                   # for curvature penalty
library(combinat)                    # for curvature penalty
library(lattice)                     # for visualization
library(Rcpp)                        # incorporate C++

source("./spline_functions.R")       # external R functions for P-splines
source("./multigrid_functions.R")    # external R functions for multigrid method
sourceCpp("./rcpp_functions.cpp")    # external C++ functions

#####--------------------------------------------
##### simple test data

set.seed(111)
x1 <- seq(0,1,length=100)
x2 <- seq(0,1,length=100)
X <- cbind(x1,x2)
P <- dim(X)[2]
X_grid <- expand.grid(x1,x2)
n <- dim(X_grid)[1]
t <- sapply( 1:n, function(i) -16*( (sum(X_grid[i,]^2) /length(X_grid[i,])) -0.5) )
fx <- 1 / ( 1 + exp(t) )
e <- rnorm(n,0,0.1)
y <- fx + e

#####--------------------------------------------
##### spline system

G <- 5                                          # number of grids for multigrid
m <- lapply( 1:G, function(g) rep(2^g-1,P) )    # number of spline knots
q <- rep(3, P)                                  # spline degree
Omega <- lapply(1:P, function(p) c( min(X[,p]),max(X[,p]) ) )   # underlying space
J <- lapply(1:G, function(g) m[[g]]+q+1)        # number of directionla basis functions
K <- prod(J[[G]])                               # total number of basis functions
tPhi_list <- lapply(1:G, function(g) lapply(1:P, function(p) t( bspline_matrix(X[,p], m[[g]][p], q[p] ,Omega[[p]]) ) ) )    # spline matrices
tPhiPhi_list <- lapply(1:G, function(g) lapply(1:P, function(p) tcrossprod(tPhi_list[[g]][[p]]) ) )
Psi_list <- lapply(1:G, function(g)  curvature_penalty(m[[g]], q, Omega) )   # survature penalty
b <- MVP_kronecker_rcpp(tPhi_list[[G]], y)       # right-hand side vector    # right-hand side vector
norm_b <- sqrt( sum(b^2) )
lambda <- 0.1                                   # weight of the regularization (manually)

#####--------------------------------------------
##### multigrid setup

Prol <- lapply( 2:G, function(g) lapply( 1:P, function(p) prolongation_matrix(J[[g-1]][p],J[[g]][p],q[p]) ) )     # prolongation matrices
Rest <- lapply( 2:G, function(g) lapply( 1:P, function(p) restriction_matrix(J[[g]][p],J[[g-1]][p],q[p]) ) )      # restriction matrices
w <- 0.5         # damping parameter for (damped) Jacobi smoother
nu <- c(6,3)     # number of pre- and post-smoothing iterations

#####--------------------------------------------
##### matrix-free MGCG method

tol <- 10^(-6)          # tolerance for stopping criterion
alpha <- rep(0, K)      # starting value
if( max(alpha^2)!=0 ){
  Aalpha <- MVP_spline(tPhiPhi_list[[G]], alpha) + lambda*MVP_penalty(Psi_list[[G]], alpha)
  r <- b-Aalpha
} else{
  r <- b
}
d <- r
z <- v_cycle(tPhiPhi_list, Psi_list, Rest, Prol, lambda, b, nu, w, alpha)    # apply MG as preconditioner
d <- z
rz <- as.numeric( crossprod(r,z) )
cat(0, sqrt(sum(r^2))/norm_b, "\n")
for(i in 1:K){         # loop of the MGCG iteration
  
  Ad <- MVP_spline(tPhiPhi_list[[G]], d) + lambda*MVP_penalty(Psi_list[[G]], d)
  t <- as.numeric( rz / crossprod(d, Ad) )
  alpha <- alpha+t*d
  rz_old <- rz
  r <- r-t*Ad
  cat(i, sqrt(sum(r^2))/norm_b, "\n")
  if(sqrt(sum(r^2))/norm_b <= tol){
    break
  }
  z <- v_cycle(tPhiPhi_list, Psi_list, Rest, Prol, lambda, r, nu, w)     # apply MG as preconditioner
  rz <- crossprod(r,z)
  beta <-  as.numeric( rz / rz_old )
  d <- z + beta*d
  
}


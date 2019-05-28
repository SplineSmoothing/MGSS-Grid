#####--------------------------------------------
##### preamble

rm(list=ls())

library(gaussquad)                   # for curvature penalty
library(combinat)                    # for curvature penalty
library(lattice)                     # for visualization
library(Rcpp)                        # incorporate C++

source("./spline_functions.R")       # external R functions for P-splines
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

m <- rep(36, P)                             # number of knots
q <- rep(3, P)                              # spline degree
Omega <- lapply(1:P, function(p) c( min(X[,p]),max(X[,p]) ) )   # underlying space
J <- m+q+1                                  # number of directionla basis functions
K <- prod(J)                                # total number of basis functions
tPhi_list <- lapply(1:P, function(p) t( bspline_matrix(X[,p], m[p], q[p] ,Omega[[p]]) ) )     # spline matrices
tPhiPhi_list <- lapply(1:P, function(p) tcrossprod(tPhi_list[[p]]) )
Psi_list <- curvature_penalty(m, q, Omega)  # curvature penalty
b <- MVP_kronecker_rcpp(tPhi_list, y)       # right-hand side vector
norm_b <- sqrt( sum(b^2) )
lambda <- 1                             # weight of the regularization (manually)

#####--------------------------------------------
##### matrix-free CG method

tol <- 10^(-6)                              # tolerance for stopping criterion
alpha <- rep(0, K)                          # starting value
if( max(alpha^2)!=0 ){
  Aalpha <-MVP_kronecker_rcpp(tPhiPhi_list, alpha) + lambda*rowSums( sapply( 1:length(Psi_list), function(p) MVP_kronecker_rcpp(Psi_list[[p]], alpha) ) )
  r <- b-Aalpha
} else{
  r <- b
}
d <- r
r_square <- crossprod(r)
cat(0, sqrt(r_square)/norm_b, "\n")
for(i in 1:K){                              # loop of the CG-iteration
  
  Ad <- MVP_kronecker_rcpp(tPhiPhi_list, d) + lambda*rowSums( sapply( 1:length(Psi_list), function(p) MVP_kronecker_rcpp(Psi_list[[p]], d) ) )
  t <- as.numeric( r_square / crossprod(d, Ad) )
  alpha <- alpha+t*d
  r_new <- r-t*Ad
  r_new_square <- crossprod(r_new)
  cat(i, sqrt(r_new_square)/norm_b, "\n")
  if(sqrt(r_new_square)/norm_b < tol){   
    break
  }
  beta <- as.numeric( r_new_square / r_square )   
  r <- r_new                                      
  r_square <- r_new_square                        
  d <- r+beta*d                 
  
}


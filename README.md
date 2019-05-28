# MGSS_grid: A Multigrid Spline Smoothing Toolbox for Grid Data
In [Siebenborn and Wagner, 2019](https://arxiv.org/abs/1901.00654) a matrix-free multigrid preconditioned CG-method is proposed to perform the penalized spline method in more than the usual `P<=2` dimensions.
The related implementations (cf. [github](https://github.com/SplineSmoothing/MGSS)) are based on exploiting the special structure of the underlying spline basis.
In the case of grid data, this structure drastically symplifies, which can be exploited to speed up the CG and the MGCG method.

## Manual
The manuals for the matrix-free CG-method (CG.R) and the matrix-free MGCG-method (MGCG.R), both with grid data, are provided.

### CG
After selecting the spline parameters, the transposed B-spline basis matrix and the curavture penalty were assembled for each spatial direction `p=1,...,P`
```R
tPhi_list <- lapply(1:P, function(p) t( bspline_matrix(X[,p], m[p], q[p] ,Omega[[p]]) ) )     # spline matrices
Psi_list <- curvature_penalty(m, q, Omega)                                                    # curvature penalty
b <- MVP_kronecker_rcpp(tPhi_list, y)                                                         # right-hand side vector
```
but since the data are gridded, the structure of the spline basis matrices can be incorporated as Kroncker product of via
```R 
tPhiPhi_list <- lapply(1:P, function(p) tcrossprod(tPhi_list[[p]]) )
```
The coefficients of the spline basis functions are determined via the solution of a linear system which is achieved by the CG-method, where the matrix-vector products with the coefficient matrix `A` are performed in a matrix-free manner.
The matrix `A` is given as a sum of Kroncker-product matrices wherefore the matrix-vector product is computed as follows
```R
Ad <- MVP_kronecker_rcpp(tPhiPhi_list, d) + lambda*rowSums( sapply( 1:length(Psi_list), function(p) MVP_kronecker_rcpp(Psi_list[[p]], d) ) )
```

### MGCG
Also the MGCG method profits from the simplified data structure by replacing every operation with the spline matrices by the corresponding Kroncker products.
For example
```R
z <- v_cycle(tPhiPhi_list, Psi_list, Rest, Prol, lambda, r, nu, w)
```
applies the `v_cycle` directly to the components `tPhiPhi_list` of the Kronecker-product.


## Application

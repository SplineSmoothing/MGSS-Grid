# MGSS-Grid: A Multigrid Spline Smoothing Toolbox for Grid Data
In [Siebenborn and Wagner, 2019](https://arxiv.org/abs/1901.00654) a matrix-free multigrid preconditioned CG-method is proposed to perform the penalized spline method in more than the usual `P<=2` dimensions.
The related implementations (cf. [github](https://github.com/SplineSmoothing/MGSS)) are based on exploiting the special structure of the underlying spline basis.
In the case of grid data, this structure drastically symplyfies, which can be exploited to speed up the CG and the MGCG method.

## Manual

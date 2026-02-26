# BSplineDensity.jl

This package provides a Julia implementation of the linear B-spline estimator for univariate data proposed by Kirkby et al. ([2021](#kirkby2021galerkin)).

The implementation is based on the Matlab [BSplineDensity](https://github.com/jkirkby3/BsplineDensity/) library, and is licensed under an MIT license as a result. Note that the implementation provided here only supports linear basis functions, and that the only option for selecting the bandwidth is least squares cross validation. For more options regarding bandwidth selection and B-spline orders, see the original implementation linked above.

## Installation
This package is not in the Julia General Registry. Install it directly from GitHub:
```julia
using Pkg
Pkg.add(url = "https://github.com/oskarhs/BSplineDensity.jl")
```

## Basic usage
The following code snippet illustrates the basic use of the package:
```julia
using BSplineDensity
using Distributions

# Simulate some data from the infamous "Claw" density
n = 10^6
d_true = MixtureModel(
    vcat(Normal(0, 1), [Normal(0.5*j, 0.1) for j in -2:2]),
    [0.5, 0.1, 0.1, 0.1, 0.1, 0.1]
)
x = rand(d_true, n)

# Fit a LinearBSplineDensity:
lsd = fit(LinearBSplineDensity, x)

# Compute the pdf of the fit at an arbitrary point or grid:
pdf(lsd, 0.2)
pdf(lsd, -3:0.01:3)

# We can also compute the cdf and the quantile functions:
cdf(lsd, -3:0.01:3)
quantile(lsd, 0.001:0.001:0.999)
```

A major advantage of this estimator is that it can be computed very quickly, even for large sample sizes. For instance, the call to fit in the above code snippet typically finishes in about half a second on my laptop (after compilation).

## Contributing
I am not planning to add any more functionality to this package at this point in time. If you are interested in extending the package, then feel free to submit a PR!

## References

<a name="kirkby2021galerkin"></a> Kirkby, J. L., Leitao, Á., and Nguyen, D. (2021). Nonparametric density estimation and bandwidth selection with B-spline bases: A novel Galerkin method. _Computational Statistics & Data Analysis_ **159**.
doi: [10.1016/j.csda.2021.107202](https://doi.org/10.1016/j.csda.2021.107202).
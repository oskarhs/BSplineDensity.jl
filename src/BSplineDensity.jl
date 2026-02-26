module BSplineDensity

using Distributions
using Interpolations
using LinearAlgebra
using StatsBase

include("galerkin.jl")
include("bandwidth.jl")

struct LinearBSplineDensity{I}
    etp::I
    function LinearBSplineDensity(grid, density)
        etp = linear_interpolation(grid, density; extrapolation_bc=0.0)
        return new{typeof(etp)}(etp)
    end
end

Distributions.pdf(lsd::LinearBSplineDensity, t::Union{Real, AbstractVector{<:Real}}) = lsd.etp(t)

function _cdf_grid(lsd::LinearBSplineDensity)
    grid    = lsd.etp.itp.ranges[1]
    d_grid  = step(grid)
    density = pdf(lsd, grid)

    cum_density = vcat(
        0.0,
        0.5 * d_grid * density[1],
        cumsum((density[1:end-1] .+ density[2:end]) .* (d_grid / 2))
    )
    cum_density = vcat(min.(cum_density, 1.0-eps()), 1.0)

    cum_grid = vcat(grid[1] - d_grid, grid, grid[end] + d_grid)

    return cum_grid, cum_density
end

function Distributions.cdf(lsd::LinearBSplineDensity, t::Union{Real, AbstractVector{<:Real}})
    cum_grid, cum_density = _cdf_grid(lsd)
    cdf_itp = interpolate(cum_grid, cum_density, SteffenMonotonicInterpolation())
    cdf_etp = extrapolate(cdf_itp, Flat())
    return cdf_etp(t)
end

function Distributions.quantile(lsd::LinearBSplineDensity, p::Union{Real, AbstractVector{<:Real}})
    cum_grid, cum_density = _cdf_grid(lsd)
    qt_itp = interpolate(cum_density, cum_grid, SteffenMonotonicInterpolation())
    qt_etp = extrapolate(qt_itp, Throw())
    return qt_etp(p)
end

function StatsBase.fit(::Type{LinearBSplineDensity}, x::AbstractVector{<:Real}; enforce_positive::Bool=true)
    N_min = 2^5
    bandwidth = BSplines_bandwidth(x)
    N, xmin, xmax = BSplines_params_by_h(x, bandwidth, N_min)
    grid = xmin .+ (0:N-1)*bandwidth
    
    beta = BSplines_coefficients_Galerkin(x, N, 1/bandwidth, xmin; enforce_positive=enforce_positive)
    h = grid[2] - grid[1]
    a = 1/h
    return LinearBSplineDensity(grid, sqrt(a) * beta)   
end

export LinearBSplineDensity
export fit, pdf, cdf

end # module

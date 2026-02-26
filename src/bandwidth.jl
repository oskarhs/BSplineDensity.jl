function BSplines_bandwidth(x::AbstractVector{<:Real})
    minimum_x, maximum_x = extrema(x)
    Range = maximum_x - minimum_x
    bandwidth = _optimize(1e-3*Range, 1e-1*Range) do h
        return BSplines_LSCV_function(h, x)
    end
    return bandwidth
end

function BSplines_LSCV(h::Real, beta::AbstractVector{<:Real}, x::AbstractVector{<:Real}, xmin::Real)
    N = length(beta)
    n = length(x)

    a    = 1 / h
    sqta = sqrt(a)

    # CV_h: tridiagonal Gram matrix product
    CV_h  = (1/3) * (beta[1]^2 + beta[N]^2)
    CV_h += (2/3) * sum(beta[k]^2 for k in 2:N-1)
    CV_h += (1/3) * sum(beta[k] * beta[k+1] for k in 1:N-1)

    # Bin samples
    K   = floor.(Int, (x .- xmin) .* a .+ 1)
    Lam = @. (x - xmin) * a + 1 - K

    cons  = sqta * sqrt(3) / n   # simplified: 3/sqrt(3) == sqrt(3)
    cons2 = (sqrt(3) - 2) * 2
    S_N   = 0.0

    for i in 1:n
        k = K[i]
        if 1 <= k <= N-1
            lam  = Lam[i]
            S_N += (1 - lam) * beta[k] + lam * beta[k+1] -
                   cons * (1 - 2*lam + 2*lam^2 + cons2 * lam * (1 - lam))
        end
    end
    S_N = (sqta / (n - 1)) * S_N

    return CV_h - 2 * S_N
end

function BSplines_params_by_h(x::AbstractVector{<:Real}, h::Real, N_min::Int=2^5)
    minimum_x, maximum_x = extrema(x)
    Range = maximum_x - minimum_x
    MIN = minimum_x - Range/2
    MAX = maximum_x + Range/2

    N = ceil(Int, (MAX - MIN) / h)
    N = max(N_min, N)

    d = ((N - 1) * h - (MAX - MIN)) / 2
    MAX = MAX + d
    MIN = MIN - d

    return N, MIN, MAX
end


function BSplines_LSCV_function(h::Real, x::AbstractVector{<:Real})

	a = 1/h
	N, xmin, xmax = BSplines_params_by_h(x, h)

	beta = BSplines_coefficients_Galerkin(x, N, a, xmin)

	CV_h = BSplines_LSCV(h, beta, x, xmin)

    return CV_h
end

# Goldern section search to find the optimum of LSCV(h)
function _optimize(
    f::F,
    x_lower::T,
    x_upper::T;
    max_iter::Int = 1000,
    rel_tol::Union{Real, Nothing} = 1e-4, # Matches MATLAB `fminbnd` defaults
    abs_tol::Union{Real, Nothing} = 1e-4  # Matches MATLAB `fminbnd` defaults
) where {F<:Function, T<:Real}

    x_lower < x_upper || error("x_lower must be less than x_upper")

    rtol = T(something(rel_tol, sqrt(eps(T))))
    atol = T(something(abs_tol, eps(T)))

    invphi::T = 0.5 * (sqrt(5) - 1)
    invphisq::T = 0.5 * (3 - sqrt(5))

    a, b = x_lower, x_upper
    h = b - a
    c = a + invphisq * h
    d = a + invphi   * h
    fc, fd = f(c), f(d)

    for _ in 1:max_iter
        if fc < fd
            m = (a + d) / 2
            (d - a) <= 2 * (atol + rtol * m) && return m
            b = d
            d, fd = c, fc
            h = b - a
            c = a + invphisq * h
            fc = f(c)
        else
            m = (c + b) / 2
            (b - c) <= 2 * (atol + rtol * m) && return m
            a = c
            c, fc = d, fd
            h = b - a
            d = a + invphi * h
            fd = f(d)
        end
    end

    error("Reached maximum number of iterations without convergence.")
end
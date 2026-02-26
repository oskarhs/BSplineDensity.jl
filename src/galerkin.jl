# function make_beta_positive(beta::AbstractVector{<:Real})
#     N = length(beta)
#     beta[1]   = max(eps(), beta[1])
#     beta[end] = max(eps(), beta[end])

#     for k in 2:N-1
#         if beta[k] < eps()
#             diff      = (eps() - beta[k]) / 2
#             beta[k]   = eps()
#             beta[k-1] = max(eps(), beta[k-1] - diff)
#             beta[k+1] = max(eps(), beta[k+1] - diff)
#         end
#     end

#     return beta
# end


# function BSplines_coefficients_Galerkin(
#     x::AbstractVector{<:Real},
#     N::Int,
#     a::Real,
#     xmin::Real;
#     enforce_positive::Bool=true
# )
#     G = Tridiagonal(
#         fill(1/6, N-1),                 # sub-diagonal
#         vcat(1/3, fill(2/3, N-2), 1/3), # diagonal
#         fill(1/6, N-1)                  # super-diagonal
#     )

#     theta = zeros(N + 1)  # extra element at end avoids boundary check in loop

#     K   = floor.(Int, (x .- xmin) .* a .+ 1)
#     XK  = xmin .+ (K .- 1) ./ a
#     Lam = (x .- XK) .* a

#     for k in 1:N
#         Lamk       = Lam[K .== k]
#         theta[k]   += sum(1 .- Lamk)
#         theta[k+1] += sum(Lamk)
#     end

#     n     = length(x)
#     theta = (a^0.5 / n) .* view(theta, 1:N)
#     beta  = G \ theta

#     if enforce_positive
#         beta = make_beta_positive(beta)
#         beta = beta .* (a^0.5 / sum(beta))
#     end

#     return beta
# end

function BSplines_coefficients_Galerkin(
    x::AbstractVector{<:Real},
    N::Int,
    a::Real,
    xmin::Real;
    enforce_positive::Bool=true
)
    G = Tridiagonal(
        fill(1/6, N-1),
        vcat(1/3, fill(2/3, N-2), 1/3),
        fill(1/6, N-1)
    )

    theta = zeros(N + 1)
    K     = floor.(Int, (x .- xmin) .* a .+ 1)
    Lam   = @. (x - xmin) * a + 1 - K

    for i in eachindex(x)
        k = K[i]
        if 1 <= k <= N
            theta[k]   += 1 - Lam[i]
            theta[k+1] +=     Lam[i]
        end
    end

    sqta = sqrt(a)
    invn = sqta / length(x)
    @. theta[1:N] *= invn        # scale in-place, no allocation
    beta = G \ view(theta, 1:N)  # view avoids copying theta

    if enforce_positive
        make_beta_positive!(beta)
        beta .*= sqta / sum(beta)  # in-place scale
    end

    return beta
end

function make_beta_positive!(beta::AbstractVector{<:Real}, eps::Real=1e-18)
    beta[1]   = max(eps, beta[1])
    beta[end] = max(eps, beta[end])
    for k in 2:length(beta)-1
        if beta[k] < eps
            diff      = (eps - beta[k]) / 2
            beta[k]   = eps
            beta[k-1] = max(eps, beta[k-1] - diff)
            beta[k+1] = max(eps, beta[k+1] - diff)
        end
    end
    return beta
end
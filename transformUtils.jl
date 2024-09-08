module TransformUtils
using DSP, AbstractFFTs, LinearAlgebra, LoopVectorization, Normalization

export TransformUtils


function elementWiseMult(
    series::Array{Vector{Float64}},
    n::Int
)::Vector{Float64}

    crossEspectrum::Vector{Float64} = ones(Float64, length(rfftfreq(n)) - 1)
    for s in series
        dft = rfft(s)
        crossEspectrum = crossEspectrum .* abs.(dft[2:length(dft)])
    end

    return crossEspectrum
end

# https://brianmcfee.net/dstbook-site/content/ch10-convtheorem/ConvolutionTheorem.html
function circConvt!(N::Vector{T}, K::Vector{T}, BS::Int) where {T<:Real}
    # Optimized for the case the kernel is in N (Shorter)
    lenN = length(N)
    lenN < length(K) && return circConvt!(K, N, BS)
    if lenN > BS
        error("BS must be >= the length of N")
    elseif length(K) > BS
        error("BS must be >= the length of K")
    end

    Y = Vector{T}(undef, BS)
    for ii in eachindex(Y)
        sumVal = zero(T)
        for kk ∈ eachindex(K)
            oa = (ii >= kk) ? N[ii-kk+1] : N[ii-kk+lenN]
            prev = sumVal + (K[kk] * oa)
            # sumVal = prev
            prev === Inf ? break : sumVal = prev
        end
        Y[ii] = sumVal
        N[ii] = sumVal
    end

    # for ii in eachindex(N)
    #     sumVal = zero(T)
    #     for kk ∈ eachindex(K)
    #         oa = (ii >= kk) ? N[ii-kk+1] : N[ii-kk+N]
    #         prev = sumVal + (K[kk] * oa)
    #         (prev === Nan || prev === Inf) ? break : sumVal = prev
    #     end
    #     N[ii] = sumVal
    # end
    return Y
end


function blockConv!(X::Vector{T}, H::Vector{T}) where {T<:Real}
    # Optimized for the case the kernel is in N (Shorter)
    Nx::Int = length(X)
    M::Int = length(H)
    Nx < M && return blockConv!(H, X)

    N = floor(Int, M + (M / 2))
    # number of zeros to pad
    M1::Int = M - 1
    # Block size
    L::Int = N - M1
    # Number of blocks
    K = floor(Int, (Nx + M1 - 1) / L)

    x::Vector{T} = reduce(vcat, [zeros(M1), X, zeros(N - 1)])

    Y = Matrix{T}(undef, 0, 0)

    @turbo for k in (0:K-1)
        xk = x[(k*L+1):(k*L+N)]
        Y[k, :] = circConvt!(xk, H, N)
    end
    Y = transpose(Y[:, M:N])
    return transpose(Y[:])

end


end
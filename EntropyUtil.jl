

module EntropyUtil

export EntropyUtil

function entropy(x::Vector{Float64})::Vector{Float64}
    """
    Calculates the entropy of a array.

    Parameters:
        x: The probabilities p0 and p1.

    Returns:
        entropy_values: The calculated entropy.

    """
    # non_zero_indices:Tuple[np.int64] = np.nonzero(x)
    # return -(x[non_zero_indices] * np.log2(x[non_zero_indices]))
    return -(x[x!=0] * np.log2(x[x!=0]))
end

function maxEntropy(
    kmers::Dict{String,Int}
)::Tuple{Float64,Int,Int,Vector{Float64}}

    """
    Kapur's threshold method.

    Kapur, J.N., P.K. Sahoo, and A.K.C. Wong. “A New Method for Gray-Level
    Picture Thresholding Using the Entropy of the Histogram.”
    Computer Vision, Graphics, and Image Processing 29,
    no. 3 (March 1985): 273–85. [DOI](https://doi.org/10.1016/0734-189X(85)90125-2).


    Parameters:
        kmers: frequency of occurrence for each kmers.

    Returns:
        maxEntropy: the maximum entropy value.
        threshold: the threshold value.
        frequency: the frequency at which the threshold occurs in the histogram.
        entropyCurve: all calculated entropies.

    """

    data::Vector{Int} = collect(values(kmers))

    totalPixel::Int = sum(data)

    descendingData::Vector{Int} = sort(data, rev=true)
    probs::Vector{Float64} = descendingData / totalPixel

    p0::Vector{Float64} = cumsum(probs)
    p1::Vector{Float64} = reverse(cumsum(probs))

    h0::Vector{Float64} = entropy(p0)
    h1::Vector{Float64} = entropy(p1)
    entropyCurve::Vector{Float64} = [h0_i + h1_i for (h0_i, h1_i) in zip(h0, h1)]

    maxEntropy::Float64 = maximum(entropyCurve)
    threshold::Int = argmax(entropyCurve)
    frequency::Int = descendingData[threshold-1]

    return (maxEntropy, threshold, frequency, entropyCurve)
end
end
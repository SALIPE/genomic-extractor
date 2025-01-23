

module EntropyUtil
include("KmerUtils.jl")

using .KmerUtils

export EntropyUtil

function _entropy(probs::Vector{Float64})::Float64
    """
    Calculates the Shannon entropy of a array.
    https://www.sciencedirect.com/topics/engineering/shannon-entropy 
    """
    return sum([-p * log2(p) for p in probs if p > 0.0])
end

function shannonEntropy(kmers::Dict{String,Int})::Float64
    data::Vector{Int} = collect(values(kmers))
    total = sum(data)
    probs::Vector{Float64} = [x / total for x in data]

    return _entropy(probs)
end

function maxEntropy(
    kmers::Dict{String,Int}
)::Float64

    """
    Kapur's threshold method.

    Kapur, J.N., P.K. Sahoo, and A.K.C. Wong. “A New Method for Gray-Level
    Picture Thresholding Using the Entropy of the Histogram.”
    Computer Vision, Graphics, and Image Processing 29,
    no. 3 (March 1985): 273–85. [DOI](https://doi.org/10.1016/0734-189X(85)90125-2).


    Parameters:
        kmers: frequency of occurrence for each kmers.

    Returns:
        frequency: the frequency at which the threshold occurs in the histogram.
    """
    data::Vector{Int} = collect(values(kmers))
    total = sum(data)
    sort!(data, rev=true)
    normalized_data::Vector{Float64} = [x / total for x in data]

    # Calculate the entropy curve
    entropy_curve = [
        let
            # Region A: probs[1:s]
            p_a = sum(normalized_data[1:s])
            h_a = p_a > 0.0 ?
                  _entropy([x / p_a for x in normalized_data[1:s]]) : 0.0

            # Region B: probs[s+1:end]
            p_b = sum(normalized_data[s+1:end])
            h_b = p_b > 0.0 ?
                  _entropy([x / p_b for x in normalized_data[s+1:end]]) : 0.0

            h_a + h_b
        end for s in 1:(length(normalized_data)-1)
    ]

    # Find the index of the maximum entropy and the corresponding frequency
    # that index will represent the treshold value
    max_entropy_idx::Int = argmax(entropy_curve)
    frequency::Int = data[max_entropy_idx]

    return entropy_curve[max_entropy_idx]
end

#Mount entropy sginal by window slide
function mountEntropyByWndw(
    wndwSize::Int,
    sequence::String)::Vector{Float64}

    index = 1
    step::Int8 = 1
    seqlen = length(sequence)

    entropy_points::Int = seqlen - wndwSize + 1
    entropyX::Vector{Float64} = zeros(Float64, entropy_points)

    pos::Int = 1
    while (index + wndwSize - 1) <= seqlen
        windown = sequence[index:index+wndwSize-1]
        kmers::Dict{String,Int} = KmerUtils.cCountKmers(windown)
        entropyX[pos] = shannonEntropy(kmers)

        index += step
        pos += 1
    end
    return entropyX
end
end

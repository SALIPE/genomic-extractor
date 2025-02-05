module KmerUtils
include("DataIO.jl")

using StatsBase,
    .DataIO

export KmerUtils



function cCountKmers(
    sequence::String,
    word::Int=3,
    step::Int=1)::Dict{String,Int}

    kmerDict = Dict{String,Int}()
    index = 1

    @inbounds while (index + word - 1) <= length(sequence)
        key = sequence[index:index+word-1]
        if haskey(kmerDict, key)
            kmerDict[key] += 1
        else
            kmerDict[key] = 1
        end
        index += step
    end

    return kmerDict
end


function selectKmers(
    threshold::Int,
    kmers::Dict{String,Int})::Dict{String,Int}

    for (k, v) in kmers
        if v > threshold
            selected_kmers[k] = 1
        end
    end

    return selected_kmers
end


function kmers_difference(
    seqKmers::Dict{String,Int},
    refKmers::Dict{String,Int}
)::Vector{String}

    seq = countmap(keys(seqKmers))
    ref = countmap(keys(refKmers))

    return collect(keys(filter((k, v) -> v > 0, seq .- ref)))
end


function kmersIntersections(
    seqKmers::Dict{String,Int},
    refKmers::Dict{String,Int}
)::Vector{String}

    seq = countmap(keys(seqKmers))
    ref = countmap(keys(refKmers))

    return collect(intersect(keys(seq), keys(ref)))
end
end

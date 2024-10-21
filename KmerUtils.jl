module KmerUtils
include("DataIO.jl")


using StatsBase,
    .DataIO

export KmerUtils

function kmersFrequencies(
    sequence::String,
    word::Int,
    step::Int)::Dict{String,Int}
    mapping = Dict('A' => 0, 'T' => 1, 'C' => 2, 'G' => 3,
        'a' => 0, 't' => 1, 'c' => 2, 'g' => 3)

    kmers = Dict{String,Int}()
    # Function to map characters to integers (A=0, T=1, C=2, G=3)

    numeric_sequence = [get(mapping, char, 4) for char in sequence]
    index = 1
    while (index + word - 1) <= length(numeric_sequence)
        DataIO.progressBar!(index + word - 1, length(numeric_sequence))

        # Extracting numeric substring (can be vectorized)
        key_numeric = numeric_sequence[index:index+word-1]

        if all(x -> x in 0:3, key_numeric)
            # Convert numeric key back to string for dictionary lookup
            key = String([['A', 'T', 'C', 'G'][n+1] for n in key_numeric])

            if haskey(kmers, key)
                kmers[key] += 1
            else
                kmers[key] = 1
            end
        end

        index += step
    end

    return kmers
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

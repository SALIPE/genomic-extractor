module KmerUtils
include("dataIO.jl")
using .DataIO
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
end

module KmerUtils
include("DataIO.jl")


using StatsBase,
    .DataIO

export KmerUtils

kmer_lock = ReentrantLock()

function _countKmers!(
    kmerDict::Dict{String,Int},
    sequence::Vector{Int},
    word::Int,
    step::Int)

    index = 1
    seqlen = length(sequence)

    while (index + word - 1) <= seqlen
        # DataIO.progressBar!(index + word - 1, seqlen)
        # Extracting numeric substring (can be vectorized)
        key_numeric = sequence[index:index+word-1]
        if all(x -> x in 0:3, key_numeric)
            # Convert numeric key back to string for dictionary lookup
            key = String([['A', 'T', 'C', 'G'][n+1] for n in key_numeric])

            if haskey(kmerDict, key)
                Threads.@lock kmer_lock kmerDict[key] += 1
            else
                Threads.@lock kmer_lock kmerDict[key] = 1
            end
        end

        index += step
    end

end

function kmersFrequencies(
    sequence::String,
    word::Int,
    step::Int)::Dict{String,Int}
    mapping = Dict('A' => 0, 'T' => 1, 'C' => 2, 'G' => 3,
        'a' => 0, 't' => 1, 'c' => 2, 'g' => 3)

    threadnum = Threads.nthreads()
    println(threadnum)

    kmers = Dict{String,Int}()
    # Function to map characters to integers (A=0, T=1, C=2, G=3)

    numeric_sequence::Vector{Int} = [get(mapping, char, 4) for char in sequence]
    seqlen = length(numeric_sequence)

    if threadnum == 1
        _countKmers!(kmers, numeric_sequence, word, step)
    else

        chunknum::Int = threadnum
        chunkSize::Int = div(seqlen, chunknum)
        tasks = Vector(undef, chunknum)

        for i in range(1, chunknum)
            if i == 1
                tasks[i] = Threads.@spawn _countKmers!(kmers, numeric_sequence[1:chunkSize], word, step)
            elseif i == chunknum
                tasks[i] = Threads.@spawn _countKmers!(kmers, numeric_sequence[(seqlen-chunkSize)-word+step:end], word, step)
            else
                tasks[i] = Threads.@spawn _countKmers!(kmers, numeric_sequence[(chunkSize*(i-1))+1-word+step:chunkSize*i], word, step)
            end
        end

        fetch.(tasks)
    end
    # count1 = Threads.@spawn _countKmers!(kmers, numeric_sequence[1:half], word, step)
    # count2 = Threads.@spawn _countKmers!(kmers, numeric_sequence[half-word-step:seqlen], word, step)

    # fetch(count1)
    # fetch(count2)

    # while (index + word - 1) <= seqlen
    #     # Extracting numeric substring (can be vectorized)
    #     key_numeric = numeric_sequence[index:index+word-1]

    #     if all(x -> x in 0:3, key_numeric)
    #         # Convert numeric key back to string for dictionary lookup
    #         key = String([['A', 'T', 'C', 'G'][n+1] for n in key_numeric])

    #         if haskey(kmers, key)
    #             kmers[key] += 1
    #         else
    #             kmers[key] = 1
    #         end
    #     end

    #     index += step
    # end

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

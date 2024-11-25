module KmerUtils
include("DataIO.jl")

using StatsBase,
    .DataIO

export KmerUtils

kmer_lock = ReentrantLock()

# For parallel execution - need improve to reuse the code below
function _countKmers!(
    kmerDict::Dict{String,Int},
    sequence::Vector{Int},
    word::Int,
    step::Int)

    index = 1
    seqlen = length(sequence)

    @inbounds while (index + word - 1) <= seqlen
        key_numeric = sequence[index:index+word-1]
        if all(x -> x in 0:3, key_numeric)
            # Convert numeric key back to string for dictionary lookup
            key = String([['A', 'T', 'C', 'G'][n+1] for n in key_numeric])
            lock(kmer_lock) do
                if haskey(kmerDict, key)
                    kmerDict[key] += 1
                else
                    kmerDict[key] = 1
                end
            end
        end
        index += step
    end
end

function _countKmers(
    sequence::Vector{Int},
    word::Int,
    step::Int)::Dict{String,Int}

    kmerDict = Dict{String,Int}()
    index = 1
    seqlen = length(sequence)

    @inbounds while (index + word - 1) <= seqlen
        key_numeric = sequence[index:index+word-1]
        if all(x -> x in 0:3, key_numeric)
            key = String([['A', 'T', 'C', 'G', 'N'][n+1] for n in key_numeric])
            if haskey(kmerDict, key)
                kmerDict[key] += 1
            else
                kmerDict[key] = 1
            end
        end
        index += step
    end

    return kmerDict
end

function kmersFrequencies(
    numeric_sequence::Vector{Int},
    word::Int=3,
    step::Int=1,
    parallel::Bool=false,
    numthreads::Int=1)::Dict{String,Int}

    threadnum = parallel ? Threads.nthreads() : numthreads

    if threadnum == 1
        return _countKmers(numeric_sequence, word, step)
    else
        kmers = Dict{String,Int}()
        seqlen = length(numeric_sequence)
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
        return kmers
    end

end

function cCountKmers(
    sequence::String,
    word::Int=3,
    step::Int=1)::Dict{String,Int}

    kmerDict = Dict{String,Int}()
    index = 1
    seqlen = length(sequence)

    @inbounds while (index + word - 1) <= seqlen
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

function kmersFrequencies(
    sequence::String,
    word::Int,
    step::Int)::Dict{String,Int}

    mapping = Dict('A' => 0, 'T' => 1, 'C' => 2, 'G' => 3,
        'a' => 0, 't' => 1, 'c' => 2, 'g' => 3)
    # Function to map characters to integers (A=0, T=1, C=2, G=3)
    numeric_sequence::Vector{Int} = [get(mapping, char, 4) for char in sequence]
    return kmersFrequencies(numeric_sequence, word, step)
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

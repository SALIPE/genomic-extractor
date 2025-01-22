module KmerUtils
include("DataIO.jl")

using StatsBase,
    .DataIO

export KmerUtils

kmer_lock = ReentrantLock()

function _countKmers!(
    kmerDict::Dict{String,Int},
    sequence::String,
    word::Int,
    step::Int)

    index = 1

    @inbounds while (index + word - 1) <= length(sequence)
        key = sequence[index:index+word-1]
        lock(kmer_lock) do
            if haskey(kmerDict, key)
                kmerDict[key] += 1
            else
                kmerDict[key] = 1
            end
        end
        index += step
    end
end



function kmersFrequencies(
    sequence::String,
    word::Int=3,
    step::Int=1,
    parallel::Bool=false,
    numthreads::Int=1)::Dict{String,Int}

    threadnum = parallel ? Threads.nthreads() : numthreads

    if threadnum == 1
        return cCountKmers(sequence, word, step)
    else
        kmers = Dict{String,Int}()
        seqlen = length(sequence)
        chunknum::Int = threadnum
        chunkSize::Int = div(seqlen, chunknum)
        tasks = Vector(undef, chunknum)

        for i in range(1, chunknum)
            if i == 1
                tasks[i] = Threads.@spawn _countKmers!(kmers, sequence[1:chunkSize], word, step)
            elseif i == chunknum
                tasks[i] = Threads.@spawn _countKmers!(kmers, sequence[(seqlen-chunkSize)-word+step:end], word, step)
            else
                tasks[i] = Threads.@spawn _countKmers!(kmers, sequence[(chunkSize*(i-1))+1-word+step:chunkSize*i], word, step)
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

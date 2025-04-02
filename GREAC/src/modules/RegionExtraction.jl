module RegionExtraction

include("DataIO.jl")
include("TransformUtils.jl")

using .DataIO,
    .TransformUtils,
    Serialization,
    FASTX,
    AbstractFFTs,
    FFTW,
    FLoops

export RegionExtraction


function regionsConjuction(
    variantDirPath::String,
    wnwPercent::Float32,
    groupName::String
)::Vector{Tuple{Int,Int}}

    variantDirs::Vector{String} = readdir(variantDirPath)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    hit_region::Union{Nothing,BitArray} = nothing

    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache("$cachdir/$(variant)_outmask.dat")

        marked = cache[2][2]

        if (isnothing(hit_region))
            hit_region = marked
        else
            for i in eachindex(marked)
                if (ismissing(hit_region[i]))
                    hit_region[i] = marked[i]
                elseif (marked[i] && !hit_region[i])
                    hit_region[i] = true
                end

            end
        end
    end

    extracted_regions = Vector{Tuple{Int,Int}}()

    init_pos = 1
    current::Bool = false

    for i in eachindex(hit_region)
        if (hit_region[i] && !current)
            current = true
            init_pos = i
        elseif (!hit_region[i] && current)
            push!(extracted_regions, (init_pos, i - 1))
            init_pos = i
            current = false
        end
    end

    if (current)
        push!(extracted_regions, (init_pos, length(hit_region)))
    end

    return extracted_regions

end


# Function to extract discriminatives features from each class
function extractFeaturesTemplate(
    wnwPercent::Float32,
    groupName::String,
    outputDir::Union{Nothing,String},
    variantDirPath::String,
    histogramThreshold::Float16=Float16(0.5))
    @info Threads.nthreads()
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    try
        mkpath(cachdir)
    catch e
        @error "create cache direcotry failed" exception = (e, catch_backtrace())
    end

    variantDirs::Vector{String} = readdir(variantDirPath)


    outputs = Vector{Tuple{String,Tuple{Vector{UInt16},BitArray}}}(undef, length(variantDirs))
    varKmer = Dict{String,Vector{String}}()

    @simd for variant in variantDirs
        varKmer[variant] = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
    end

    exclusiveKmers::Dict{String,Vector{String}} = findExclusiveElements(varKmer)


    @inbounds for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Processing $variant")
        cache_path = "$cachdir/$(variant)_outmask.dat"

        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray}}} = DataIO.load_cache(cache_path)

        if !isnothing(cache)
            @info "Using cached data from $cache_path"
            outputs[v] = cache
        else

            sequences::Vector{String} = DataIO.loadStringSequences("$variantDirPath/$variant/$variant.fasta")
            minSeqLength::UInt64 = minimum(map(length, sequences))
            wnwSize::UInt64 = ceil(UInt64, minSeqLength * wnwPercent)

            data::Tuple{String,Tuple{Vector{UInt16},BitArray}} = (variant, _wndwExlcusiveKmersHistogram(exclusiveKmers[variant], wnwSize, sequences, histogramThreshold))
            outputs[v] = data
            DataIO.save_cache(cache_path, data)
        end
    end

end


function findExclusiveElements(class_dict::Dict{String,Vector{String}})::Dict{String,Vector{String}}
    # Track which classes each string appears in
    presence = Dict{String,Set{String}}()

    # First pass: Record class presence for each unique string
    for (class_name, strings) in class_dict
        unique_strings = unique(strings)
        for str in unique_strings
            if haskey(presence, str)
                push!(presence[str], class_name)
            else
                presence[str] = Set([class_name])
            end
        end
    end

    # Second pass: Find exclusive elements for each class
    exclusive_dict = Dict{String,Vector{String}}()

    for (class_name, strings) in class_dict
        unique_elements = unique(strings)
        exclusive = filter(str -> presence[str] == Set([class_name]), unique_elements)
        exclusive_dict[class_name] = exclusive
    end

    return exclusive_dict
end

#=
    Feature extraction based on exclusive kmers.
    Extraction process:
    k-mers list -> run windown slide -> count histogram for most freq -> extract regions

=#
function wndwExlcusiveKmersHistogram(
    exclusiveKmers::Vector{String},
    wndwSize::UInt64,
    sequences::Vector{String},
    histogramThreshold::Float16
)::Tuple{Vector{UInt16},BitArray}

    @assert all(≤(wndwSize), length.(exclusiveKmers)) "All k-mers must be ≤ window size"

    patterns = [Base.Fix1(occursinKmer, kmer) for kmer in exclusiveKmers]

    byte_seqs = [codeunits(s) for s in sequences]
    maxSeqLen = maximum(length, sequences)
    total_windows = maxSeqLen - wndwSize + 1

    # num_threads = Threads.nthreads()
    # basesize = max(1, length(sequences) ÷ num_threads)

    @floop for seq in byte_seqs
        seq_len = length(seq)
        seq_windows = seq_len - wndwSize + 1
        seq_hist = zeros(UInt16, seq_windows)

        for initPos in 1:seq_windows
            endPos = initPos + wndwSize - 1
            for pattern in patterns
                if pattern(seq[initPos:endPos])
                    seq_hist[initPos] += 1
                    break
                end
            end
        end

        padded_hist = zeros(UInt16, total_windows)
        valid_range = 1:length(seq_hist)
        padded_hist[valid_range] = seq_hist

        @reduce(
            histogram = zeros(UInt16, total_windows) .+ padded_hist
        )
    end

    marked = falses(maxSeqLen)
    for i in eachindex(histogram)
        if histogram[i] / length(byte_seqs) > histogramThreshold
            marked[i:i+wndwSize-1] .= true
        end
    end
    return histogram, marked
end

function _wndwExlcusiveKmersHistogram(
    exclusiveKmers::Vector{String},
    wndwSize::UInt64,
    sequences::Vector{String},
    histogramThreshold::Float16
)::Tuple{Vector{UInt16},BitArray}

    @assert all(≤(wndwSize), length.(exclusiveKmers)) "All k-mers must be ≤ window size"

    k_len = length(exclusiveKmers[1])
    byte_seqs = [codeunits(s) for s in sequences]
    maxSeqLen = maximum(length, sequences)
    total_windows = maxSeqLen - wndwSize + 1

    @floop for seq in byte_seqs
        padded_hist = zeros(UInt16, total_windows)

        positions = getOccursin(String(seq), exclusiveKmers)

        for i in positions

            iniPos::Int = i - (wndwSize - k_len)
            iniPos = iniPos > 1 ? iniPos : 1
            endPos = i > total_windows ? total_windows : i
            padded_hist[iniPos:endPos] = ones(UInt16, endPos - iniPos + 1)
        end

        @reduce(
            histogram = zeros(UInt16, total_windows) .+ padded_hist
        )
    end

    marked = falses(maxSeqLen)
    for i in eachindex(histogram)
        if histogram[i] / length(byte_seqs) > histogramThreshold
            marked[i:i+wndwSize-1] .= true
        end
    end
    return histogram, marked
end


function wndwSequencesKmersHistogram(
    kmerset::Set{String},
    wndwSize::UInt16,
    sequences::Vector{String},
)::Vector{Vector{UInt64}}

    kmer_lengths = length.(kmerset)
    @assert all(≤(wndwSize), kmer_lengths) "All k-mers must be ≤ window size"

    patterns = [Base.Fix1(occursinKmer, kmer) for kmer in kmerset]

    byte_seqs = [codeunits(s) for s in sequences]
    maxSeqLen = maximum(length, sequences)
    total_windows = maxSeqLen - wndwSize + 1

    x_signals = Vector{Vector{UInt64}}(undef, length(sequences))

    @floop for (i, seq) in enumerate(byte_seqs)
        seq_windows = length(seq) - wndwSize + 1
        seq_hist = zeros(UInt64, seq_windows)

        for initPos in 1:seq_windows
            endPos = initPos + wndwSize - 1
            for pattern in patterns
                if pattern(seq[initPos:endPos])
                    seq_hist[initPos] += 1
                end
            end
        end

        padded_hist = zeros(UInt64, total_windows)
        padded_hist[1:length(seq_hist)] = seq_hist

        x_signals[i] = padded_hist
    end

    return x_signals

end


#=
Count kmers appearance in the region segment, doing a codeunits regex
=#
function countPatterns(
    seqWindow::SubArray,
    kmers::Vector{String})::UInt16

    patterns = [Base.Fix1(occursinKmer, kmer) for kmer in kmers]
    count::UInt16 = 0

    @floop for pattern in patterns
        if pattern(seqWindow)
            @reduce count += 1
        end
    end
    return count
end

function getOccursin(
    sequence::String,
    kmerset::Vector{String})::Vector{Int}

    occurrences_pos = Set{Int}()
    @inbounds for i in eachindex(kmerset)
        regex = Regex("\\Q$(kmerset[i])\\E")
        match_positions = [m.offset for m in eachmatch(regex, sequence)]
        union!(occurrences_pos, Set(match_positions))
    end
    return collect(occurrences_pos)
end

function occursinKmer(
    kmer::String,
    windowBuffer::Union{SubArray,Vector{UInt8}}
)::Bool
    return occursin(Regex("\\Q$kmer\\E"), String(windowBuffer))
end

#=
Codeunits "regex", beacuse is a byte comparison 
=#
function occursinKmerBit(
    kmer::Base.CodeUnits,
    windowBuffer::Union{SubArray,Vector{UInt8}}
)::Bool
    wlen = length(windowBuffer)
    klen = length(kmer)

    @inbounds for i in 1:(wlen-klen+1)
        match = true
        for j in 1:klen
            (windowBuffer[i+j-1] ≠ kmer[j]) && (match = false; break)
        end
        match && return true
    end
    return false

end

end



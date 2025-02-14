module Model

include("DataIO.jl")
include("RRMUtils.jl")
include("Classification.jl")

using .DataIO,
    .RRM,
    .Classification,
    FLoops,
    FASTX

export Model

function trainModel(
    wnwPercent::Float32,
    outputDir::String,
    variantDirPath::String)

    cachdir::String = "$(pwd())/.project_cache"

    try
        mkdir(cachdir)
    catch e
        @error "create cache direcotry failed" exception = (e, catch_backtrace())
    end

    variantDirs::Vector{String} = readdir(variantDirPath)
    @show Threads.nthreads()

    outputs = Vector{Tuple{String,Tuple{Vector{UInt16},BitArray},Vector{String}}}(undef, length(variantDirs))
    varKmer = Dict{String,Vector{String}}()

    @simd for variant in variantDirs
        varKmer[variant] = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
    end

    exclusiveKmers::Dict{String,Vector{String}} = findExclusiveElements(varKmer)


    for v in eachindex(variantDirs)
        variant::String = variantDirs[v]
        println("Processing $variant")
        cache_path = "$cachdir/$(variant)_outmask.dat"

        cache::Union{Nothing,Tuple{String,Tuple{Vector{UInt16},BitArray},Vector{String}}} = DataIO.load_cache(cache_path)

        if !isnothing(cache)
            @info "Using cached data from $cache_path"
            outputs[v] = cache
        else

            sequences = String[]
            for record in open(FASTAReader, "$variantDirPath/$variant/$variant.fasta")
                push!(sequences, sequence(String, record))
            end
            minSeqLength::UInt16 = minimum(map(length, sequences))
            wnwSize::UInt16 = ceil(UInt16, minSeqLength * wnwPercent)

            data::Tuple{String,Tuple{Vector{UInt16},BitArray},Vector{String}} = (variant, wndwExlcusiveKmersHistogram(exclusiveKmers[variant], wnwSize, sequences), sequences)
            outputs[v] = data
            DataIO.save_cache(cache_path, data)
        end
        println("Finish Processing $variant")
    end

    #Structure defined as {Variant Name , (Regions Marked, Fourier Coefficients, Exclusive Kmers)}
    # trainedModel = Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}()
    trainedModel = Dict([(variant, (BitArray{1}(), Vector{Vector{Float64}}(), String[])) for variant in variantDirs])

    for (variant, (_, marked), sequences) in outputs
        fourierCoefficients = Vector{Vector{Float64}}()

        minSeqLength::UInt16 = minimum(map(length, sequences))
        limitedMark::BitArray = marked[1:minSeqLength]
        start = 0
        current = false

        for (i, bit) in enumerate(limitedMark)
            if bit && !current
                start = i
                current = true
            elseif !bit && current
                cross = RRM.getFourierCoefficient(
                    [str[start:i-1] for str in sequences],
                    i - 1 - start
                )
                push!(fourierCoefficients, cross)
                current = false
            end
        end
        if current
            cross = RRM.getFourierCoefficient(
                [str[start:minSeqLength] for str in sequences],
                minSeqLength - start
            )
            push!(fourierCoefficients, cross)
        end
        trainedModel[variant] = (marked, fourierCoefficients, exclusiveKmers[variant])
    end

    DataIO.save_cache("$cachdir/trained_model.dat", trainedModel)

    for (variant, (histogram, marked)) in outputs
        plt = plot(histogram, title="Exclusive Kmers Histogram - $wnwPercent", dpi=300)
        png(plt, "$outputDir/$variant")
        plt = plot(marked, title="Exclusive Kmers Marked - $wnwPercent", dpi=300)
        png(plt, "$outputDir/$(variant)_reg")
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

#extract from k-mers
# k-mers list -> run windown slide -> count histogram for most freq -> extract regions
function wndwExlcusiveKmersHistogram(
    exclusiveKmers::Vector{String},
    wndwSize::UInt16,
    sequences::Vector{String},
)::Tuple{Vector{UInt16},BitArray}

    kmer_lengths = length.(exclusiveKmers)
    @assert all(≤(wndwSize), kmer_lengths) "All k-mers must be ≤ window size"

    patterns = [Base.Fix1(Classification.occursinKmerBit, codeunits(kmer)) for kmer in exclusiveKmers]

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
                    # break
                end
            end
        end

        padded_hist = zeros(UInt16, total_windows)
        valid_range = 1:length(seq_hist)
        padded_hist[valid_range] = seq_hist

        @reduce(
            histogram = ones(UInt16, total_windows) .* padded_hist,
            # histogram = zeros(UInt16, total_windows) .+ padded_hist,
        )
    end

    marked = falses(maxSeqLen)

    for i in eachindex(histogram)
        if histogram[i] > 0
            marked[i:i+wndwSize-1] .= true
        end
    end
    return histogram, marked

end


end



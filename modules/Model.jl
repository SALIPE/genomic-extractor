module Model

include("DataIO.jl")
include("TransformUtils.jl")
include("Classification.jl")

using .DataIO,
    .TransformUtils,
    .Classification,
    DecisionTree,
    ScikitLearn,
    AbstractFFTs,
    FLoops

export Model

mutable struct ClassificationBlockRegressor
    class::String
    models::Vector{Tuple{Int,Int,Vector{Float64},Any}}
end

mutable struct ModelClassBlockStruct
    blockmodelchain::Array{ClassificationBlockRegressor}
end

function createModelBlock(
    blockClassName::String,
    regionsFourierCoefs::Vector{Tuple{Int,Int,Vector{Float64}}},
    trueSequences::AbstractArray,
    falseSequences::AbstractArray,
    nTree
)::ClassificationBlockRegressor

    block_model = ClassificationBlockRegressor(blockClassName, Vector{Tuple{Int,Int,Any}}())

    # Pre-allocate reusable containers
    X_buffer = Vector{Float64}[]  # Reused for each region
    labels_buffer = Float16[]     # Reused for each region

    for (initIdx, endIdx, cross) in regionsFourierCoefs
        freIndexes = findall(x -> x > 0, cross)
        isempty(freIndexes) && continue

        # Fix critical min/max frequency calculation
        minFreq = minimum(freIndexes)
        maxFreq = maximum(freIndexes)
        region_length = endIdx - initIdx + 1


        # Precompute FFT plan for this region
        fft_plan = plan_rfft(Vector{Float64}(undef, region_length))

        # Reuse buffers to reduce allocations
        empty!(X_buffer)
        empty!(labels_buffer)
        sizehint!(X_buffer, length(falseSequences) + length(trueSequences))
        sizehint!(labels_buffer, length(falseSequences) + length(trueSequences))

        # Process sequences using view-based access
        process_sequences!(X_buffer, labels_buffer, falseSequences, initIdx, endIdx, fft_plan, minFreq, maxFreq, 0.0f16)
        process_sequences!(X_buffer, labels_buffer, trueSequences, initIdx, endIdx, fft_plan, minFreq, maxFreq, 1.0f16)

        # isempty(X_buffer) && continue

        # Convert to matrix without intermediate allocations
        X_matrix = Matrix{Float64}(undef, length(X_buffer), length(X_buffer[1]))
        for i in eachindex(X_buffer)
            X_matrix[i, :] = X_buffer[i]
        end

        regionBlockModel = RandomForestRegressor(n_trees=nTree)
        DecisionTree.fit!(regionBlockModel, X_matrix, labels_buffer)

        push!(block_model.models, (initIdx, endIdx, cross, regionBlockModel))
    end

    return block_model
end

function process_sequences!(X, labels, sequences, initIdx, endIdx, fft_plan, minFreq, maxFreq, label)
    @inbounds for seq in sequences
        lastindex(seq) < endIdx && continue
        # Use views to avoid substring allocations
        subseq = @view seq[initIdx:endIdx]
        num = DataIO.sequence2AminNumSerie(subseq)

        # Compute FFT using pre-planned transform
        dft = abs.(fft_plan * num)
        dft_trimmed = @view dft[2:end]  # Skip DC component


        push!(X, dft_trimmed[minFreq:maxFreq])
        push!(labels, label)
    end
end

function trainModel(
    classes::Dict{String,Vector{String}},
    model::Dict{String,Tuple{BitArray,Vector{Tuple{Int,Int,Vector{Float64}}},Vector{String}}},
    nTrees=20
)::ModelClassBlockStruct

    modelBlckStruct = ModelClassBlockStruct(Vector{ClassificationBlockRegressor}(undef, length(classes)))
    class_names = collect(keys(classes))

    # Precompute all false sequences once
    total_sequences = sum(length(v) for v in values(classes))
    class_index = Dict{String,Vector{Int}}()
    all_seqs = Vector{String}(undef, total_sequences)

    # Build global sequence index
    idx = 1
    for (cls, seqs) in classes
        class_index[cls] = idx:(idx+length(seqs)-1)
        all_seqs[idx:(idx+length(seqs)-1)] = seqs
        idx += length(seqs)
    end

    @inbounds for (i, cls) in enumerate(class_names)
        true_range = class_index[cls]
        false_ranges = [v for (k, v) in class_index if k != cls]

        # Use view-based selection for memory efficiency
        trueSequences = @view all_seqs[true_range]
        falseSequences = vcat([@view all_seqs[r] for r in false_ranges]...)

        modelBlckStruct.blockmodelchain[i] = createModelBlock(
            cls,
            model[cls][2],
            trueSequences,
            falseSequences,
            nTrees
        )
    end

    return modelBlckStruct
end


# Functio to extract discriminatives features from each class
function extractFeaturesTemplate(
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


    @inbounds for v in eachindex(variantDirs)
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
    trainedModel = Dict([(variant, (BitArray{1}(), Vector{Tuple{Int,Int,Vector{Float64}}}(), String[])) for variant in variantDirs])


    # the idea here its to create a way to extract a way how the class should be using RRM, this way
    # we well have the windows and the frequence range and it magnitude (0-1) to compare.

    for (variant, (_, marked), sequences) in outputs
        fourierCoefficients = Vector{Tuple{Int,Int,Vector{Float64}}}()

        minSeqLength::UInt16 = minimum(map(length, sequences))
        limitedMark::BitArray = marked[1:minSeqLength]
        start = 0
        current = false

        for (i, bit) in enumerate(limitedMark)
            if bit && !current
                start = i
                current = true
            elseif !bit && current
                cross = TransformUtils.getFourierCoefficient([str[start:i-1] for str in sequences])
                push!(fourierCoefficients, (start, i - 1, cross))
                current = false
            end
        end
        if current
            cross = TransformUtils.getFourierCoefficient([str[start:minSeqLength] for str in sequences])
            push!(fourierCoefficients, (start, minSeqLength, cross))
        end
        trainedModel[variant] = (marked, fourierCoefficients, exclusiveKmers[variant])
    end

    DataIO.save_cache("$cachdir/trained_model.dat", trainedModel)

    # for (variant, (histogram, marked)) in outputs
    #     plt = plot(histogram, title="Exclusive Kmers Histogram - $wnwPercent", dpi=300)
    #     png(plt, "$outputDir/$variant")
    #     plt = plot(marked, title="Exclusive Kmers Marked - $wnwPercent", dpi=300)
    #     png(plt, "$outputDir/$(variant)_reg")
    # end

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
            histogram = zeros(UInt16, total_windows) .+ padded_hist
        )
    end

    marked = falses(maxSeqLen)
    # regionCount = count(i -> i > 0, histogram)
    # probabilities = Vector{Float64}(undef, regionCount)

    # p_idx = one(Int)
    for i in eachindex(histogram)
        if histogram[i] > 0
            marked[i:i+wndwSize-1] .= true
            # probabilities[p_idx] = countin[i] / length(sequences)
            # p_idx += 1
        end
    end
    return histogram, marked

end


end



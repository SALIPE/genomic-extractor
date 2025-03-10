module Model

include("DataIO.jl")
include("TransformUtils.jl")

using .DataIO,
    .TransformUtils,
    Serialization,
    DecisionTree,
    FASTX,
    AbstractFFTs,
    FFTW,
    FLoops

export Model

mutable struct ClassificationBlockRegressor
    class::String
    models::Vector{Tuple{Int,Int,Vector{Float64},Any}}
end

mutable struct ModelClassBlockStruct
    blockmodelchain::Array{ClassificationBlockRegressor}
end

#=
Create blocks for region classification for each class and mask pre-extracted data structure
=#
function createModelBlock(
    blockClassName::String,
    regionsFourierCoefs::Vector{Tuple{Int,Int,Vector{Float64}}},
    trueSequences::AbstractVector{String},
    falseSequences::AbstractVector{String},
    nTree
)::ClassificationBlockRegressor

    block_model = ClassificationBlockRegressor(blockClassName, Vector{Tuple{Int,Int,Any}}())


    for (initIdx, endIdx, cross) in regionsFourierCoefs
        freIndexes = findall(>(0), cross)
        isempty(freIndexes) && continue

        minFreq::Int = minimum(freIndexes)
        maxFreq::Int = maximum(freIndexes)

        X = Vector{Vector{float(Float32)}}()
        labels = Vector{Float16}()

        #extract all the fft from each region
        for seq in falseSequences
            lastindex(seq) < endIdx && continue

            subseq = @view seq[initIdx:endIdx]
            num = DataIO.sequence2AminNumSerie(subseq)
            dft = abs.(rfft(num)[2:end])

            push!(X, dft[minFreq:maxFreq])
            push!(labels, zero(Float16))

        end
        for seq in trueSequences
            lastindex(seq) < endIdx && continue

            subseq = @view seq[initIdx:endIdx]
            num = DataIO.sequence2AminNumSerie(subseq)
            dft = abs.(rfft(num)[2:end])

            push!(X, dft[minFreq:maxFreq])
            push!(labels, one(Float16))
        end


        X_matrix = reduce(hcat, X) |> permutedims

        regionBlockModel = RandomForestRegressor()
        fit!(regionBlockModel, X_matrix, labels)
        push!(block_model.models, (initIdx, endIdx, cross, regionBlockModel))
    end
    return block_model
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
    outputDir::Union{Nothing,String},
    variantDirPath::String)

    cachdir::String = "$(pwd())/.project_cache/$wnwPercent"

    try
        mkpath(cachdir)
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
        #     plt = plot(histogram, title="Exclusive Kmers Histogram - $wnwPercent", dpi=300)
        #     png(plt, "$outputDir/$variant")
        #     plt = plot(marked, title="Exclusive Kmers Marked - $wnwPercent", dpi=300)
        #     png(plt, "$outputDir/$(variant)_reg")

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

    DataIO.save_cache("$cachdir/extracted_features.dat", trainedModel)

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
    wndwSize::UInt16,
    sequences::Vector{String},
)::Tuple{Vector{UInt16},BitArray}

    kmer_lengths = length.(exclusiveKmers)
    @assert all(≤(wndwSize), kmer_lengths) "All k-mers must be ≤ window size"

    patterns = [Base.Fix1(occursinKmerBit, codeunits(kmer)) for kmer in exclusiveKmers]

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


function wndwSequencesKmersHistogram(
    kmerset::Set{String},
    wndwSize::UInt16,
    sequences::Vector{String},
)::Vector{Vector{UInt64}}

    kmer_lengths = length.(kmerset)
    @assert all(≤(wndwSize), kmer_lengths) "All k-mers must be ≤ window size"

    patterns = [Base.Fix1(occursinKmerBit, codeunits(kmer)) for kmer in kmerset]

    byte_seqs = [codeunits(s) for s in sequences]
    maxSeqLen = maximum(length, sequences)
    total_windows = maxSeqLen - wndwSize + 1

    # kmer_total::Int = length(kmerset)
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
        # x_signals[i] = Float64.(padded_hist) ./ kmer_total
    end

    return x_signals

end

function classifyInput(
    inputSequence::AbstractString,
    wnwPercent::Float32,
    outputfilename::Union{Nothing,AbstractString}=nothing
)

    @show outputfilename
    modelCachedFile = "$(pwd())/.project_cache/$wnwPercent/extracted_features.dat"

    model::Union{Nothing,Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}} = DataIO.load_cache(modelCachedFile)

    if !isnothing(model)
        @info "Using model from cached data from $modelCachedFile"
        classifyInput(inputSequence, model, outputfilename)
    else
        error("Model not found in cached files!")
    end


end

#=
Classification function based on the counters of exclusive kmers appearance.

model structure description (this is the extracted features): Dict{VariantName, (Marked, Fourier Coefficients, Kmers)}
=#
function classifyInput(
    inputSequence::Base.CodeUnits,
    model::Dict{String,Tuple{BitArray,Vector{Tuple{Int,Int,Vector{Float64}}},Vector{String}}},
    outputfilename::Union{Nothing,String}
)

    report = Dict{String,Vector{Tuple{Tuple{UInt16,UInt16},UInt16}}}()
    for (key, (marked, _, kmers)) in model
        inputlen = minimum(length, [inputSequence, marked])
        report[key] = Vector{Tuple{Tuple{UInt16,UInt16},UInt16}}()
        limitedMark::BitArray = marked[1:inputlen]
        start = 0
        current = false

        for (i, bit) in enumerate(limitedMark)
            if bit && !current
                start = i
                current = true
            elseif !bit && current
                current = false
                count = @views countPatterns(inputSequence[start:i-1], kmers)

                push!(report[key], ((start, i - 1), count))
            end
        end
        if current
            count = @views countPatterns(inputSequence[start:inputlen], kmers)
            push!(report[key], ((start, inputlen), count))
        end
    end

    if !isnothing(outputfilename)
        open(outputfilename, "w") do file
            for (var, regions) in report
                write(file, "\n\n########### $(uppercase(var)) ############")
                write(file, "\nTotal Exclusive Kmers: $(length(model[var][3]))")
                write(file, "\n####################################\n")
                for ((initPos, endPos), count) in regions
                    write(file, "\nWindow Position( $initPos - $endPos ): \nKmer Count: $count")
                end
            end
        end
    end

    classifications = Dict{String,Float16}()

    for (var, regions) in report
        predictions::BitArray = []
        for ((initPos, endPos), count) in regions
            push!(predictions, count > 0)
        end
        classifications[var] = count(i -> i, predictions) / length(predictions)
    end

    maxPercent = maximum(x -> x[2], classifications)
    return findfirst(x -> x == maxPercent, classifications), maxPercent, classifications

end


function classifySequence(
    blockModel::ModelClassBlockStruct,
    inputSequence::String)

    classifications = Dict{String,Float16}()
    for block in blockModel.blockmodelchain

        predictions::BitArray = []
        for (initIdx, endIdx, cross, innerModel) in block.models
            lastindex(inputSequence) < endIdx && continue

            freIndexes = findall(>(0), cross)
            isempty(freIndexes) && continue

            minFreq::Int = minimum(freIndexes)
            maxFreq::Int = maximum(freIndexes)

            subseq = @view inputSequence[initIdx:endIdx]
            num = DataIO.sequence2AminNumSerie(subseq)

            dft = abs.(rfft(num)[2:end])
            norm = dft[minFreq:maxFreq]

            input = reduce(hcat, [norm]) |> permutedims
            raw_output = predict(innerModel, input)[1]

            push!(predictions, raw_output >= 0.5f0)
        end
        if length(predictions) == 0
            @error "Non discriminative points for classification on class $(block.class)"
        else
            classifications[block.class] = count(i -> i, predictions) / length(predictions)
        end

    end
    maxPercent = maximum(x -> x[2], classifications)
    classification::Union{String,Nothing} = findfirst(x -> x == maxPercent, classifications)

    if isnothing(classification)
        error("Insufficient data for classification")
    end
    return classification, maxPercent, classifications
end


#=
Count kmers appearance in the region segment, doing a codeunits regex
=#
function countPatterns(
    seqWindow::SubArray,
    kmers::Vector{String})::UInt16

    patterns = [Base.Fix1(occursinKmerBit, codeunits(kmer)) for kmer in kmers]
    count::UInt16 = 0

    @floop for pattern in patterns
        if pattern(seqWindow)
            @reduce count += 1
        end
    end
    return count
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



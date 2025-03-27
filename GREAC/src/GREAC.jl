
module GREAC

include("modules/DataIO.jl")
include("modules/EntropyUtil.jl")
include("modules/TransformUtils.jl")
include("modules/ConvergenceAnalysis.jl")
include("modules/RegionExtraction.jl")
include("modules/ClassificationModel.jl")

using FLoops,
    FASTX,
    Plots,
    LinearAlgebra,
    Normalization,
    Statistics,
    ArgParse,
    .DataIO,
    .RegionExtraction,
    .ClassificationModel,
    .TransformUtils,
    .EntropyUtil,
    .ConvergenceAnalysis


export FastaSplitter

function findPosWndw(
    positions::Vector{Int},
    wndwSize::Int,
    entropyX::Vector{Float64})::Vector{Float64}

    regions = zeros(Float64, length(entropyX))

    for pos in positions
        initIdx::Int = pos + 1 - wndwSize
        initIdx = initIdx <= 0 ? 1 : initIdx

        endIdx::Int = initIdx + wndwSize - 1
        endIdx = endIdx > length(entropyX) ? length(entropyX) : endIdx
        regions[initIdx:endIdx] .= 1
    end
    # como cada ponto representa uma região e queremos ver as regiões que conteem aquela posição
    # precisamos analisae da seguinte maneira: um ponto em X = pos->pos+wndwSize - 1
    # portanto precisamos selecionar os indices que tem o ponto, do que terminam com aquele ponto
    # ate os que começam naquele ponto, vulgo a intersecção das regiões que tem a posição
    return regions
end

function histogramPosWndw(
    positions::Vector{Int16},
    wndwSize::Int16,
    signalLen::Int32)::Vector{Int16}

    posHistogram = zeros(Int16, signalLen)

    for pos in positions
        initIdx::Int16 = pos + 1 - wndwSize
        initIdx = initIdx <= 0 ? 1 : initIdx

        # endIdx::Int16 = initIdx + wndwSize - 1
        endIdx::Int16 = pos > signalLen ? signalLen : pos
        posHistogram[initIdx:endIdx] .+= 1
    end
    return posHistogram
end

# Exec the window process for one FASTA file
function validateEntropyWindow(
    positions::Vector{Int16},
    wnwPercent::Float16,
    variantName::String,
    output::String,
    variantFilePath::String,
    validatePos::Bool=true
)

    sequences::Array{String} = []
    for record in open(FASTAReader, variantFilePath)
        seq::String = sequence(String, record)
        push!(sequences, seq)
    end
    # plt = plot(title="$variantName Regions - $(wnwPercent*100)%", dpi=300)

    entropy_signals = Vector{Vector{Float64}}(undef, length(sequences))
    # Processo para encontrar valores de entropia por região do genoma
    for (s, seqs) in enumerate(sequences)
        slideWndw::Int = ceil(Int, length(seqs) * wnwPercent)
        y::Vector{Float64} = EntropyUtil.mountEntropyByWndw(slideWndw, seqs)
        entropy_signals[s] = y

        # plot!(x, y, label="Entropy-value", xlims=lim)
    end
    # distances = ConvergenceAnalysis.euclidean_distance(entropy_signals)
    # ylen = length(distances)
    # x = range(1, ylen)
    # lim = [0, ylen]
    # if validatePos
    #     regions::Vector{Int} = histogramPosWndw(positions, ceil(Int, ylen * wnwPercent), ylen)
    #     plot!(twinx(), x, regions,
    #         label="Frequency",
    #         seriestype=:bar,
    #         linecolor=nothing,
    #         xlims=lim)
    # end
    # plot!(x, distances, label="Distance-value", xlims=lim)
    # png(plt, output)
end


function computeEntropySignal(
    filePath::String,
    wnwPercent::Float32,
)::Vector{Vector{Float64}}

    sequences::Array{String} = []

    for record in open(FASTAReader, filePath)
        seq::String = sequence(String, record)
        push!(sequences, seq)
    end

    num_sequences = length(sequences)
    entropy_signals = Vector{Vector{Float64}}(undef, num_sequences)

    if num_sequences > 1
        @floop for s = 1:num_sequences
            seq = sequences[s]
            y = EntropyUtil.mountEntropyByWndw(
                ceil(Int, length(seq) * wnwPercent), seq)
            @inbounds entropy_signals[s] = y
        end
        return entropy_signals
    else
        seq = sequences[1]
        y = EntropyUtil.mountEntropyByWndw(
            ceil(Int, length(seq) * wnwPercent), seq)
        entropy_signals[1] = y
        return entropy_signals
    end

end

#=
     Do a comparison betweeen variant class regions based on the euclidian distance
     the comparison should be a explanation to why use a consensus signal for class classfication evaluation
 =#
function compareVariantClassPerDistance(
    wnwPercent::Float32,
    output::String,
    variantDirPath::String,
    positions::Vector{Int16}=Int16[],
    valPositions::Bool=false
)

    files::Vector{String} = readdir(variantDirPath)

    consensusSignals = Vector{Tuple{String,Vector{Float64}}}(undef, length(files))

    @show Threads.nthreads()

    @floop ThreadedEx() for (i, file) in enumerate(files)
        println("Processing $file")

        entropy_signals = computeEntropySignal("$variantDirPath/$file", wnwPercent)

        # Signifca que é um arquivo consensus (pelo menos deveria ser)
        if length(entropy_signals) == 1
            consensusSignals[i] = (file, entropy_signals[1])
        else
            distances::Vector{Float64} = ConvergenceAnalysis.euclidean_distance(entropy_signals)
            consensusSignals[i] = (file, distances)
        end

        println("Finish Processing $file")

    end

    # plt = plot(title="Variant Classes Comparison per window size - $wnwPercent", dpi=300)

    # for (fastaFile, distances) in consensusSignals
    #     plot!(range(1, length(distances)), distances, label="$fastaFile: Distance-value")
    # end

    #png(plt, "$output/euclidian_consensus")

    filteredEntropy = TransformUtils.RRMEntropySignal(consensusSignals)
    # plt = plot(title="Signals Filtered using RRM - $wnwPercent", dpi=300)

    ylen::Int32 = minimum(map(x -> length(x[2]), filteredEntropy))
    x = range(1, ylen)
    lim = [0, ylen]

    # for (class, f) in filteredEntropy
    #     plot!(x, f[1:ylen], label=class, xlims=lim)
    # end

    regions::Vector{Int16} = Vector{Int16}()

    # if valPositions

    #     regions = histogramPosWndw(positions, ceil(Int16, ylen * wnwPercent), ylen)
    #     plot!(twinx(), x, regions,
    #         label="Frequency",
    #         seriestype=:bar,
    #         linecolor=nothing,
    #         xlims=lim)

    # end
    # png(plt, "$output/iffts")

    _, _, norm = findPeaksBetweenClasses(map(x -> x[2], consensusSignals))


    # plt = plot(title="Signal distances between points classes - $wnwPercent", dpi=300)
    # plot!(norm,
    #     linecolor=:red,
    #     xlims=lim)
    # plot!(twinx(), x, regions,
    #     label="Frequency",
    #     seriestype=:bar,
    #     linecolor=nothing,
    #     xlims=lim)

    # png(plt, "$output/distances")

end


#Função para achar quais são os picos de distancia entre as classes, a aprtir do consensus de distancia gerado 
function findPeaksBetweenClasses(signals::Vector{Vector{Float64}})
    n_classes = length(signals)

    # Tamanho de cada vetor
    numPoints = minimum(map(length, signals))

    # Matriz para armazenar distâncias ponto a ponto
    matrixDistances = zeros(Float64, numPoints, n_classes, n_classes)

    # Calcular distâncias ponto a ponto para todas as combinações de classes
    for i in 1:n_classes, j in i+1:n_classes
        for k in 1:numPoints
            matrixDistances[k, i, j] = sqrt((signals[i][k] - signals[j][k])^2)
            matrixDistances[k, j, i] = matrixDistances[k, i, j]  # Simetria
        end
    end

    # Matriz média de distâncias para cada ponto
    matrixMeamDistances = mean(matrixDistances, dims=(2, 3))[:]

    N = MinMax(matrixMeamDistances)
    norm = N(matrixMeamDistances)

    # Mediana das médias pra definir picos
    peakThreashold = maximum(matrixMeamDistances)
    # Identificar os picos de maior valor
    peaksIdx = findall(x -> x == peakThreashold, matrixMeamDistances)

    return (matrixMeamDistances, peaksIdx, norm)
end



function greacClassification(
    folderPath::String,
    outputdir::Union{Nothing,String},
    wnwPercent::Float32,
    groupName::String,
    metric::Union{Nothing,String}
)

    modelCachedFile = "$(homedir())/.project_cache/$groupName/$wnwPercent/kmers_distribution.dat"
    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache(modelCachedFile)

    confMatrix = Dict{String,Tuple{Int,Int}}()
    classification_probs = Dict{String,Vector{Tuple{String,Dict{String,Float64}}}}()

    classify = Base.Fix1(ClassificationModel.predict_raw, (model, metric))

    y_true = Vector{String}()
    y_pred = Vector{String}()

    for class in model.classes
        classeqs::Vector{Base.CodeUnits} = DataIO.loadCodeUnitsSequences("$folderPath/$class.fasta")

        @info "Classyfing $class sequences:"
        classifications = Vector{Tuple{String,Dict{String,Float64}}}(undef, length(classeqs))
        for (i, seq) in enumerate(classeqs)
            push!(y_true, class)
            input::Vector{Base.CodeUnits} = [seq]
            get_appearences = Base.Fix1(ClassificationModel.def_kmer_classes_probs, (model.regions, input))

            @floop for kmer in collect(model.kmerset)
                kmer_seq_histogram = get_appearences(kmer)

                @reduce(
                    kmer_distribution = zeros(UInt64, length(model.regions)) .+ kmer_seq_histogram
                )
            end
            seq_distribution = kmer_distribution ./ length(model.kmerset)
            classifications[i] = classify(seq_distribution)

            push!(y_pred, classifications[i][1])
        end

        classification_probs[class] = classifications
        confMatrix[class] = (count(x -> x[1] == class, classifications), length(classifications))

    end
    @info confMatrix

    results = compute_variant_metrics(model.classes, y_true, y_pred)

    println("######### Results for :$wnwPercent ###########")
    # Access results:
    println("Confusion Matrix:")
    display(results[:confusion_matrix])

    println("\nPer-Class Metrics:")
    for (variant, metrics) in results[:per_class]
        println("$variant: ", metrics)
    end

    println("\nMacro Averages: ", results[:macro])
    println("Micro Averages: ", results[:micro])
    println("Weighted Averages: ", results[:weighted])

end

function compute_variant_metrics(
    variants::Vector{String},
    y_true::Vector{String},
    y_pred::Vector{String})
    n_classes = length(variants)

    # Validate input
    length(y_true) == length(y_pred) || error("Input vectors must have the same length")
    all(v in variants for v in y_true) || error("Invalid variant in y_true")
    all(v in variants for v in y_pred) || error("Invalid variant in y_pred")

    # Create confusion matrix
    class_idx = Dict(v => i for (i, v) in enumerate(variants))
    cm = zeros(Int, n_classes, n_classes)

    for (t, p) in zip(y_true, y_pred)
        i = class_idx[t]
        j = class_idx[p]
        cm[i, j] += 1
    end

    # Calculate per-class metrics
    metrics = Dict{String,Dict}()
    total_samples = sum(cm)

    for (i, variant) in enumerate(variants)
        tp = cm[i, i]
        fp = sum(cm[:, i]) - tp
        fn = sum(cm[i, :]) - tp
        tn = total_samples - (tp + fp + fn)

        precision = (tp + fp) == 0 ? 0.0 : tp / (tp + fp)
        recall = (tp + fn) == 0 ? 0.0 : tp / (tp + fn)
        f1 = (precision + recall) ≈ 0.0 ? 0.0 : 2 * (precision * recall) / (precision + recall)
        support = sum(cm[i, :])

        metrics[variant] = Dict(
            :precision => round(precision, digits=4),
            :recall => round(recall, digits=4),
            :f1 => round(f1, digits=4),
            :support => support
        )
    end

    # Aggregate metrics
    macro_precision = mean([m[:precision] for m in values(metrics)])
    macro_recall = mean([m[:recall] for m in values(metrics)])
    macro_f1 = mean([m[:f1] for m in values(metrics)])

    micro_tp = sum(diag(cm))
    micro_precision = micro_tp / total_samples
    micro_recall = micro_precision  # Same as accuracy
    micro_f1 = micro_precision

    weighted_precision = sum([m[:precision] * m[:support] for m in values(metrics)]) / total_samples
    weighted_recall = sum([m[:recall] * m[:support] for m in values(metrics)]) / total_samples
    weighted_f1 = sum([m[:f1] * m[:support] for m in values(metrics)]) / total_samples

    return Dict(
        :confusion_matrix => cm,
        :classes => variants,
        :per_class => metrics,
        :macro => Dict(
            :precision => round(macro_precision, digits=4),
            :recall => round(macro_recall, digits=4),
            :f1 => round(macro_f1, digits=4)
        ),
        :micro => Dict(
            :precision => round(micro_precision, digits=4),
            :recall => round(micro_recall, digits=4),
            :f1 => round(micro_f1, digits=4)
        ),
        :weighted => Dict(
            :precision => round(weighted_precision, digits=4),
            :recall => round(weighted_recall, digits=4),
            :f1 => round(weighted_f1, digits=4)
        )
    )
end

function getKmersDistributinPerClass(
    wnwPercent::Float32,
    groupName::String,
    variantDirPath::String
)
    cachdir::String = "$(homedir())/.project_cache/$groupName/$wnwPercent"

    try
        mkpath(cachdir)
        @info "Cache dir created"
    catch e
        @error "create cache directory failed" exception = (e, catch_backtrace())
    end

    model::Union{Nothing,ClassificationModel.MultiClassModel} = DataIO.load_cache("$cachdir/kmers_distribution.dat")

    if !isnothing(model)
        @info "Model already processed from cached data from $cachdir"
    else

        variantDirs::Vector{String} = readdir(variantDirPath)

        kmerset = Set{String}()
        kmers_dist = Dict{String,Int}()

        @simd for variant in variantDirs
            variantKmers = DataIO.read_pickle_data("$variantDirPath/$variant/$(variant)_ExclusiveKmers.sav")
            kmers_dist[variant] = length(variantKmers)
            union!(kmerset, Set(variantKmers))
        end

        meta_data = Dict{String,Int}()
        byte_seqs = Dict{String,Vector{Base.CodeUnits}}()
        wnw_size = one(Int)
        max_seq_len = one(Int)

        for variant in variantDirs
            byte_seqs[variant] = DataIO.loadCodeUnitsSequences("$variantDirPath/$variant/$variant.fasta")

            minSeqLength::Int = minimum(map(length, byte_seqs[variant]))
            maxSeqLength::Int = maximum(map(length, byte_seqs[variant]))

            class_wnw_size = ceil(Int, minSeqLength * wnwPercent)

            if (wnw_size == one(Int) || wnw_size > class_wnw_size)
                wnw_size = class_wnw_size
            end

            if (max_seq_len == one(Int) || maxSeqLength > max_seq_len)
                max_seq_len = maxSeqLength
            end


            meta_data[variant] = length(byte_seqs[variant])
        end

        @info "Window size value: $wnw_size"

        max_seq_windows = max_seq_len - wnw_size + 1

        @info meta_data

        distribution::ClassificationModel.MultiClassModel = ClassificationModel.fitMulticlass(
            kmerset,
            kmers_dist,
            meta_data,
            byte_seqs,
            wnw_size,
            max_seq_windows,
            RegionExtraction.regionsConjuction(variantDirPath, wnwPercent, groupName))

        DataIO.save_cache("$cachdir/kmers_distribution.dat", distribution)
    end
end


function main()
    settings = ArgParseSettings(
        description="Genome Regions Extractor and Classifier",
        version="0.1",
        add_version=true,
        prog="greac"
    )

    # Global options (apply to all commands)
    @add_arg_table! settings begin
        "--no-cache"
        help = "Remove cached files"
        action = :store_true
        "--group-name"
        help = "Process Group Name"
        required = true
        arg_type = String

    end

    # Create subcommand structure
    @add_arg_table! settings begin
        ("convergence-analysis", action=:command,
            help="Perform convergence analysis between sequences")
        ("extract-model", action=:command,
            help="Extract discriminative features for classification")
        ("classify", action=:command,
            help="Classify sequences using a pre-trained model")
        ("benchmark", action=:command,
            help="Benchmark extract features model and classify creating and print confusion matrix")
        ("region-validation", action=:command,
            help="Validate genomic regions with various parameters")
    end

    # Add arguments for each subcommand
    add_convergence_analysis_args!(settings)
    add_extract_model_args!(settings)
    add_classify_args!(settings)
    add_benchmark_args!(settings)
    add_region_validation_args!(settings)

    parsed_args = parse_args(settings)


    try

        if parsed_args["no-cache"]
            rm("$(homedir())/.project_cache/$(parsed_args["group-name"])"; recursive=true, force=true)
        end

        if parsed_args["%COMMAND%"] == "convergence-analysis"
            handle_convergence_analysis(parsed_args["convergence-analysis"])
        elseif parsed_args["%COMMAND%"] == "extract-model"
            handle_extract_model(parsed_args["extract-model"])
        elseif parsed_args["%COMMAND%"] == "classify"
            handle_classify(parsed_args["classify"])
        elseif parsed_args["%COMMAND%"] == "benchmark"
            handle_benchmark(parsed_args["benchmark"], parsed_args["group-name"])
        elseif parsed_args["%COMMAND%"] == "region-validation"
            handle_region_validation(parsed_args["region-validation"])
        end
    catch e
        @error "Error processing command" exception = (e, catch_backtrace())
    end
end

function add_convergence_analysis_args!(settings)
    s = settings["convergence-analysis"]
    @add_arg_table! s begin
        "-d", "--files-directory"
        help = "Variant directory with organisms Fasta file"
        required = true
        arg_type = String
        "-w", "--window"
        help = "Slide window percent size to apply"
        arg_type = Float32
        default = 0.004
    end
end

function add_extract_model_args!(settings)
    s = settings["extract-model"]
    @add_arg_table! s begin
        "-d", "--files-directory"
        help = "Variant directory with organisms Fasta file"
        required = true
        arg_type = String
        "-w", "--window"
        help = "Slide window percent size to apply"
        arg_type = Float32
        default = 0.004
        "-o", "--output-directory"
        help = "Output directory for the model"
        arg_type = String
    end
end

function add_classify_args!(settings)
    s = settings["classify"]
    @add_arg_table! s begin
        "-w", "--window"
        help = "Sliding window percent size"
        arg_type = Float32
        required = true
        "--test-dir"
        help = "Test dataset path"
        required = true
        "-m", "--metric"
        help = "Metric used for classification"
        required = false
        "-o", "--output-directory"
        help = "Output directory for results"
        arg_type = String
    end
end

function add_benchmark_args!(settings)
    s = settings["benchmark"]
    @add_arg_table! s begin
        "-w", "--window"
        help = "Sliding window percent size"
        arg_type = Float32
        required = true
        range_tester = (x -> 0.001 < x < 0.5)
        "--train-dir"
        help = "Training dataset path"
        required = true
        "-m", "--metric"
        help = "Metric used for classification"
        required = false
        range_tester = (x -> x in ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        "--test-dir"
        help = "Test dataset path"
        required = true
    end
end

function add_region_validation_args!(settings)
    s = settings["region-validation"]
    @add_arg_table! s begin
        "-f", "--file"
        help = "Single Fasta file"
        arg_type = String
        "-d", "--files-directory"
        help = "Directory containing Fasta files"
        arg_type = String
        "--variant-name"
        help = "Variant name identifier"
        arg_type = String
        "-w", "--window"
        help = "Sliding window percent size"
        arg_type = Float32
        default = 0.004
        "-p", "--val-positions-file"
        help = "Validation positions file"
        arg_type = String
        "-o", "--output-directory"
        help = "Output directory"
        required = true
        arg_type = String
    end
end

#= Command Handlers =#

function handle_convergence_analysis(args)
    @info "Starting convergence analysis" args
    ConvergenceAnalysis.convergenceAnalysis(
        args["window"],
        args["files-directory"]
    )
end

function handle_extract_model(args)
    @info "Starting model extraction" args

    RegionExtraction.extractFeaturesTemplate(
        args["window"],
        nothing,
        args["files-directory"]
    )

    getKmersDistributinPerClass(
        args["window"],
        args["files-directory"]
    )


end

function handle_classify(args)
    @info "Starting classification"
    greacClassification(
        args["test-dir"],
        nothing,
        args["window"],
        args["group-name"],
        args["metric"]
    )
end

function handle_benchmark(args, groupName::String)
    @info "Starting benchmark" args
    @info "Starting model extraction" args
    RegionExtraction.extractFeaturesTemplate(
        args["window"],
        groupName,
        nothing,
        args["train-dir"]
    )
    getKmersDistributinPerClass(
        args["window"],
        groupName,
        args["train-dir"]
    )
    @info "Starting classification evaluation" args
    greacClassification(
        args["test-dir"],
        nothing,
        args["window"],
        groupName,
        args["metric"]
    )
end

function handle_region_validation(args)
    @info "Starting region validation" args
    positions = isnothing(args["val-positions-file"]) ? nothing :
                DataIO.readVectorFromFile(args["val-positions-file"], Int16)

    if !isnothing(args["file"])
        validateEntropyWindow(
            positions,
            args["window"],
            args["variant-name"],
            args["output-directory"],
            args["file"]
        )
    elseif !isnothing(args["files-directory"])
        compareVariantClassPerDistance(
            args["window"],
            args["output-directory"],
            args["files-directory"],
            positions,
            true
        )
    else
        error("Must specify either --file or --files-directory")
    end
end



if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end

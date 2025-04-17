module GREAC

include("modules/DataIO.jl")
include("modules/EntropyUtil.jl")
include("modules/TransformUtils.jl")
include("modules/ConvergenceAnalysis.jl")
include("modules/RegionExtraction.jl")
include("modules/ClassificationModel.jl")

using FLoops,
    FASTX,
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


export main

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

    sequences::Vector{String} = DataIO.loadStringSequences(filePath)

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

    @floop ThreadedEx() for (i, file) in enumerate(files)

        entropy_signals = computeEntropySignal("$variantDirPath/$file", wnwPercent)

        # Signifca que é um arquivo consensus (pelo menos deveria ser)
        if length(entropy_signals) == 1
            consensusSignals[i] = (file, entropy_signals[1])
        else
            distances::Vector{Float64} = ConvergenceAnalysis.euclidean_distance(entropy_signals)
            consensusSignals[i] = (file, distances)
        end
    end

    # -------------------- SIGNAL COMPARISON --------------------------------------------
    # plt = plot(title="Variant Classes Comparison per window size - $wnwPercent", dpi=300)
    # for (fastaFile, distances) in consensusSignals
    #     plot!(range(1, length(distances)), distances, label="$fastaFile: Distance-value")
    # end
    # savefig(plt, "$output/euclidian_consensus.pdf")
    # -----------------------------------------------------------------------------------


    # -------------------- SIGNAL FILTERING --------------------------------------------
    # filteredEntropy = TransformUtils.RRMEntropySignal(consensusSignals)
    # plt = plot(title="Signals Filtered using RRM - $wnwPercent", dpi=300)

    # ylen::Int32 = minimum(map(x -> length(x[2]), filteredEntropy))
    # x = range(1, ylen)
    # lim = [0, ylen]

    # for (class, f) in filteredEntropy
    #     plot!(x, f[1:ylen], label=class, xlims=lim)
    # end

    # regions::Vector{Int16} = Vector{Int16}()

    # if valPositions

    #     regions = histogramPosWndw(positions, ceil(Int16, ylen * wnwPercent), ylen)
    #     plot!(twinx(), x, regions,
    #         label="Frequency",
    #         seriestype=:bar,
    #         linecolor=nothing,
    #         xlims=lim)

    # end
    # png(plt, "$output/iffts")

    # _, _, norm = findPeaksBetweenClasses(map(x -> x[2], consensusSignals))


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
    # -----------------------------------------------------------------------------------

end


#Função para achar quais são os picos de distancia entre as classes, a partir do consensus de distancia gerado 
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

    # classification_probs = Dict{String,Vector{Tuple{String,Dict{String,Float64}}}}()

    classify = Base.Fix1(ClassificationModel.predict_raw, (model, metric))

    y_true = String[]
    y_pred = String[]
    kmerset::Vector{String} = collect(model.kmerset)

    for class in model.classes
        file_path::String = "$folderPath/$class.fasta"
        total = DataIO.countSequences(file_path)

        chunk_size = 10000
        chunk_init::Int = 1

        @info "Classyfing $class $total sequences:"
        # classifications = Vector{Tuple{String,Dict{String,Float64}}}(undef, total)
        local_y_pred = Vector{String}(undef, total)
        local_y_true = Vector{String}(undef, total)
        fill!(local_y_true, class)

        while chunk_init <= total

            chunk_end = min(chunk_init + chunk_size - 1, total)
            current_chunk_size = chunk_end - chunk_init + 1

            classeqs::Vector{Base.CodeUnits} = DataIO.loadCodeUnitsSequences(file_path, chunk_init, chunk_end)

            # inner_classifications = Vector{Tuple{String,Dict{String,Float64}}}(undef, current_chunk_size)
            inner_y_pred = Vector{String}(undef, current_chunk_size)

            @floop for local_idx in 1:current_chunk_size
                seq::Base.CodeUnits = classeqs[local_idx]

                kmer_distribution = ClassificationModel.sequence_kmer_distribution(
                    model.regions, seq, kmerset
                )
                seq_distribution = kmer_distribution ./ length(model.kmerset)

                if !iszero(sum(seq_distribution))
                    cl = classify(seq_distribution)
                    # inner_classifications[local_idx] = cl
                    inner_y_pred[local_idx] = cl[1]
                end
            end

            # classifications[chunk_init:chunk_end] = inner_classifications
            local_y_pred[chunk_init:chunk_end] = inner_y_pred
            GC.gc(false)

            @info "Chunk processed $chunk_init - $chunk_end"
            chunk_init = chunk_end + 1
        end

        append!(y_pred, local_y_pred)
        append!(y_true, local_y_true)
        # classification_probs[class] = classifications

    end

    results = compute_variant_metrics(model.classes, y_true, y_pred)
    @info results
    # RESULTS_CSV = "benchmark_results_$groupName.csv"

    # open(RESULTS_CSV, "a") do io
    #     # Write header if file is empty/new
    #     if filesize(RESULTS_CSV) == 0
    #         types = join(model.classes, ",")
    #         write(io, "wndwPercent,metric,windows,window_size,max_seq_windows,kmerset," * types * ",macro_f1,macro_precision,macro_recall,micro_f1,micro_precision,micro_recall\n")
    #     end

    #     # Format data components
    #     # cm = replace(string(results[:confusion_matrix]), "\n" => " | ")
    #     perclass = join([v[:f1] for (k, v) in results[:per_class]], ",")
    #     # Create CSV line
    #     line = join([
    #             escape_string(string(wnwPercent)),
    #             escape_string(string(metric)),
    #             length(model.regions),
    #             model.wnw_size,
    #             model.max_seq_windows,
    #             length(model.kmerset),
    #             perclass,
    #             results[:macro][:f1],
    #             results[:macro][:precision],
    #             results[:macro][:recall],
    #             results[:micro][:f1],
    #             results[:micro][:precision],
    #             results[:micro][:recall]
    #         ], ",")

    #     write(io, line * "\n")
    # end

    # println("######### Results for :$wnwPercent  - $metric ###########")
    # # Access results:
    # println("Confusion Matrix:")
    # display(results[:confusion_matrix])

    # println("\nPer-Class Metrics:")
    # for (variant, metrics) in results[:per_class]
    #     println("$variant: ", metrics)
    # end

    # println("\nMacro Averages: ", results[:macro])
    # println("Micro Averages: ", results[:micro])

end

function compute_variant_metrics(
    variants::Vector{String},
    y_true::Vector{String},
    y_pred::Vector{String}
)
    n_classes = length(variants)

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

    metrics = Dict{String,Dict}()
    total_samples = sum(cm)

    total_tp = 0
    total_fp = 0
    total_fn = 0

    for (i, variant) in enumerate(variants)
        tp = cm[i, i]
        fp = sum(cm[:, i]) - tp
        fn = sum(cm[i, :]) - tp

        # For micro metrics
        total_tp += tp
        total_fp += fp
        total_fn += fn

        precision = (tp + fp) == 0 ? 0.0 : tp / (tp + fp)
        recall = (tp + fn) == 0 ? 0.0 : tp / (tp + fn)
        f1 = (precision + recall) ≈ 0.0 ? 0.0 : 2 * ((precision * recall) / (precision + recall))
        support = sum(cm[i, :])

        metrics[variant] = Dict(
            :precision => round(precision, digits=4),
            :recall => round(recall, digits=4),
            :f1 => round(f1, digits=4),
            :support => support
        )
    end

    macro_precision = mean([m[:precision] for m in values(metrics)])
    macro_recall = mean([m[:recall] for m in values(metrics)])
    macro_f1 = mean([m[:f1] for m in values(metrics)])

    micro_precision = (total_tp + total_fp) == 0 ? 0.0 : total_tp / (total_tp + total_fp)
    micro_recall = (total_tp + total_fn) == 0 ? 0.0 : total_tp / (total_tp + total_fn)
    micro_f1 = (micro_precision + micro_recall) ≈ 0.0 ? 0.0 :
               2 * ((micro_precision * micro_recall) / (micro_precision + micro_recall))

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

        #=
        This values was used to have the total number of features when working
        with the sliding window at the entire sequences, and the sequences in dataset
        have different lengths, this way we always work with a fixed feature array, 
        considering all the information independently of the input length
        =#
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
        "-w", "--window"
        help = "Sliding window percent size"
        arg_type = Float32
        required = true
        range_tester = (x -> 0.001 < x < 0.5)

    end

    # Create subcommand structure
    @add_arg_table! settings begin
        ("convergence-analysis", action=:command,
            help="Perform convergence analysis between sequences")
        ("benchmark", action=:command,
            help="Benchmark extract features model and classify creating and print confusion matrix")
        ("region-validation", action=:command,
            help="Validate genomic regions with various parameters")
    end

    # Add arguments for each subcommand
    add_convergence_analysis_args!(settings)
    add_benchmark_args!(settings)
    add_region_validation_args!(settings)

    parsed_args = parse_args(settings)


    try

        if parsed_args["no-cache"]
            rm("$(homedir())/.project_cache/$(parsed_args["group-name"])/$(parsed_args["window"])"; recursive=true, force=true)
        end

        if parsed_args["%COMMAND%"] == "convergence-analysis"
            handle_convergence_analysis(parsed_args["convergence-analysis"])
        elseif parsed_args["%COMMAND%"] == "benchmark"
            handle_benchmark(parsed_args["benchmark"], parsed_args["group-name"], parsed_args["window"])
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

    end
end


function add_benchmark_args!(settings)
    s = settings["benchmark"]
    @add_arg_table! s begin
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

function handle_benchmark(args,
    groupName::String,
    window::Float32)
    @info "Starting benchmark" args window groupName
    @info "Starting model extraction"
    RegionExtraction.extractFeaturesTemplate(
        window,
        groupName,
        nothing,
        args["train-dir"]
    )
    getKmersDistributinPerClass(
        window,
        groupName,
        args["train-dir"]
    )
    @info "Starting classification evaluation"
    greacClassification(
        args["test-dir"],
        nothing,
        window,
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

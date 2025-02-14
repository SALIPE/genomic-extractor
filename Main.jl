
using Distributed

# Add worker processes first
addprocs(5)

@everywhere begin
    using Pkg
    Pkg.activate(".")
    Pkg.resolve()
    Pkg.instantiate()

end

@everywhere include("modules/DataIO.jl")
@everywhere include("modules/KmerUtils.jl")
@everywhere include("modules/EntropyUtil.jl")
@everywhere include("modules/TransformUtils.jl")
@everywhere include("modules/RRMUtils.jl")
@everywhere include("modules/ConvergenceAnalysis.jl")
@everywhere include("modules/Classification.jl")
@everywhere include("modules/Model.jl")

@everywhere using FLoops,
    FASTX,
    LinearAlgebra,
    Normalization,
    Statistics,
    ArgParse

@everywhere using
    .DataIO,
    .Model,
    .RRM,
    .Classification,
    .KmerUtils,
    .TransformUtils,
    .EntropyUtil,
    .ConvergenceAnalysis

begin

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
        distances = ConvergenceAnalysis.euclidean_distance(entropy_signals)
        ylen = length(distances)
        x = range(1, ylen)
        lim = [0, ylen]
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


    # Do a comparison betweeen variant class regions based on the euclidian distance
    # the comparison should be a explanation to why use a consensus signal for class classfication evaluation
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

        filteredEntropy = RRM.RRMEntropySignal(consensusSignals)
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



    function extractKmerFeaturesByFrequence(
        wnwPercent::Float32,
        outputDir::String,
        variantDirPath::String
    )

        Model.trainModel(wnwPercent, outputDir, variantDirPath)

    end

    function sequencesClassification(
        fasta::String,
        modelPath::Union{Nothing,String},
        outputdir::String
    )

        @show fasta
        sequences = Vector{Tuple{String,Base.CodeUnits}}()
        for record in open(FASTAReader, fasta)
            seq::String = sequence(String, record)
            id::String = identifier(record)
            push!(sequences, (replace(id, r"\/|\|" => "_"), codeunits(seq)))
            break
        end

        modelCachedFile = isnothing(modelPath) ? "$(pwd())/.project_cache/trained_model.dat" : modelPath

        model::Union{Nothing,Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}} = DataIO.load_cache(modelCachedFile)

        if !isnothing(model)
            @info "Using model from cached data from $modelCachedFile"

            @floop for (id, seq) in sequences
                @show id
                Classification.classifyInput(seq, model, "$outputdir/$id")
            end

        else
            error("Model not found in cached files!")
        end


    end

    function validate_region_main()
        setting = ArgParseSettings()
        @add_arg_table! setting begin
            "-f", "--file"
            help = "Single Fasta file"
            "-d", "--files-directory"
            help = "Variant directory with organisms Fasta file"
            "--variant-name"
            help = "Variant name"
            arg_type = String
            "-w", "--window"
            help = "Slide window percent size to apply"
            arg_type = Float16
            default = 0.004
            "-p", "--val-positions-file"
            help = "Positions file for validation"
            arg_type = String
            "--convergence-analysis"
            help = "Positions file for validation"
            action = :store_true
            "--extract-model"
            help = "Mount K-mers positions"
            action = :store_true
            "--classify"
            help = "Create sequence report classification"
            "-o", "--output-directory"
            help = "Output directory"
            arg_type = String
            required = true
        end

        parsed_args = parse_args(ARGS, setting)

        filePath::Union{Nothing,String} = parsed_args["file"]
        dirPath = parsed_args["files-directory"]
        windowSize::Float32 = parsed_args["window"]
        positionsBedFile = parsed_args["val-positions-file"]
        outputDirectory::String = parsed_args["output-directory"]
        varname = parsed_args["variant-name"]
        execConvAnalysis = parsed_args["convergence-analysis"]
        kmersFreq = parsed_args["extract-model"]
        classify::Union{Nothing,String} = parsed_args["classify"]

        println("Parsed args:")
        for (arg, val) in parsed_args
            println("  $arg  =>  $val")
        end

        if execConvAnalysis
            # ConvergenceAnalysis.convergenceAnalysis(windowSize, dirPath)
            #ConvergenceAnalysis.convergenceAnalysisClasses(windowSize, dirPath)
            return 0
        elseif kmersFreq
            validateKmerFrequencies(windowSize, outputDirectory, dirPath)
            return 0
        elseif !isnothing(classify)
            sequencesClassification(filePath, classify, outputDirectory)
            return 0

        elseif !isnothing(positionsBedFile)
            positions::Vector{Int16} = DataIO.readVectorFromFile(positionsBedFile, Int16)

            if !isnothing(filePath)
                # Recebe so o arquivo fasta de uma unic variant
                validateEntropyWindow(positions, windowSize, varname, outputDirectory, filePath)
            elseif !isnothing(dirPath)
                # Rece um diretório como arquivos fastas dividido por classes, cada arquivo contem N organismos
                compareVariantClassPerDistance(windowSize, outputDirectory, dirPath, positions, true)
            else
                println("Non mode selected")
            end
        else
            if !isnothing(dirPath)
                # Rece um diretório como arquivos fastas dividido por classes, cada arquivo contem N organismos
                compareVariantClassPerDistance(windowSize, outputDirectory, dirPath)
            else
                println("Non mode selected")
            end
        end
    end

    validate_region_main()

end

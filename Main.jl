
using Distributed, Pkg
Pkg.instantiate()
Pkg.activate(".")

include("DataIO.jl")
include("KmerUtils.jl")
include("EntropyUtil.jl")
include("TransformUtils.jl")
include("RRMUtils.jl")

# addprocs(4)

# @everywhere

using AbstractFFTs,
    FASTX,
    Plots,
    LoopVectorization,
    Normalization,
    ArgParse,
    .DataIO,
    .KmerUtils,
    .TransformUtils,
    .EntropyUtil

using .RRM: runRRMMethodology!
begin


    function rrm_main()
        setting = ArgParseSettings()
        @add_arg_table! setting begin
            "-w", "--window"
            help = "Slide window size to apply method in the chromossome. Values from [256, 512, 1024, 2048]"
            default = 2048
            "-r", "--range-tolerance"
            help = "Tolarence range to create ranges ex. -t 2 => [(1-15), (18-20)] turns into [(1-20)]"
            default = 100
            "-t", "--threshold"
            help = "Threshold value for the filtered signal to select regions"
            default = 0.2
            "-b", "--bed-file"
            help = "BED file with regions to analyse"
        end

        parsed_args = parse_args(ARGS, setting)

        sldWindow::Int = parsed_args["window"]
        tolerance::Int = parsed_args["range-tolerance"]
        threshold::Float64 = parsed_args["threshold"]
        regionsbedfile = parsed_args["bed-file"]

        valid_wndw = [256, 512, 1024, 2048]
        if !(sldWindow in valid_wndw)
            println("Error: Invalid window size. Choose from: ", join(valid_wndw, ", "))
            exit(1)
        end

        println("Selected slide window size: $sldWindow")
        println("Selected tolerance value: $tolerance")
        println("Selected threshold value: $threshold")

        filePath::String = "/home/salipe/Desktop/GitHub/datasets/gen_dron_car.fasta"
        dirPath::String = "/home/salipe/Desktop/GitHub/datasets/dron_iber_test"

        if regionsbedfile === nothing
            # runRRMMethodology!(filePath, sldWindow, tolerance, threshold)
        else
            regions = DataIO.readRegionsFromBed(regionsbedfile)
            runRRMMethodology!(dirPath, regions, tolerance, threshold)
        end
    end

    function max_entropy_main()
        setting = ArgParseSettings()
        @add_arg_table! setting begin
            "-f", "--file"
            help = "Fasta file"
            required = true
        end

        parsed_args = parse_args(ARGS, setting)

        fastaFile::String = parsed_args["file"]
        sequences::Array{String} = []
        for record in open(FASTAReader, fastaFile)
            seq::String = sequence(String, record)
            push!(sequences, seq)
        end
        kmers = KmerUtils.kmersFrequencies(sequences[1], 6, 1)
        @show kmers
        frequency = EntropyUtil.maxEntropy(kmers)
        @show frequency

    end

    function mountEntropyByWndw(
        wndwSize::Int,
        step::Int8,
        sequence::String)::Vector{Float64}

        index = 1

        seqlen = length(sequence)
        pos::Int16 = 1
        entropyX::Vector{Float64} = zeros(Float64, seqlen - wndwSize + 1)

        while (index + wndwSize - 1) <= seqlen
            windown = sequence[index:index+wndwSize-1]
            kmers::Dict{String,Int} = KmerUtils.cCountKmers(windown)
            entropyValue = EntropyUtil.shannonEntropy(kmers)
            entropyX[pos] = entropyValue
            pos += 1
            index += step
        end

        return entropyX
    end

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
        positions::Vector{Int},
        wndwSize::Int,
        signalLen::Int)::Vector{Int}

        posHistogram = zeros(Int, signalLen)

        for pos in positions
            initIdx::Int = pos + 1 - wndwSize
            initIdx = initIdx <= 0 ? 1 : initIdx

            endIdx::Int = initIdx + wndwSize - 1
            endIdx = endIdx > signalLen ? signalLen : endIdx
            posHistogram[initIdx:endIdx] .+= 1
        end
        return posHistogram
    end

    function validateEntropyWindow(
        positions::Vector{Int},
        wnwPercent::Float16,
        variantName::String,
        testLbl::String,
        variantFilePath::String
    )

        sequences::Array{String} = []
        for record in open(FASTAReader, variantFilePath)
            seq::String = sequence(String, record)
            push!(sequences, seq)
        end
        wndwStep::Int8 = 1
        plt = plot(title="$variantName Regions - $(wnwPercent*100)%")

        # Processo para encontrar valores de entropia por região do genoma
        for seqs in sequences
            slideWndw::Int = ceil(Int, length(seqs) * wnwPercent)
            y::Vector{Float64} = mountEntropyByWndw(slideWndw, wndwStep, seqs)

            ylen = length(y)
            x = range(1, ylen)
            lim = [0, ylen]
            regions::Vector{Int} = histogramPosWndw(positions, slideWndw, ylen)

            plot!(twinx(), x, regions,
                label="Frequency",
                seriestype=:bar,
                linecolor=nothing,
                xlims=lim)
            plot!(x, y, label="Entropy-value", xlims=lim)
        end

        png(plt, testLbl)
    end


    function validate_region_main()
        setting = ArgParseSettings()
        @add_arg_table! setting begin
            "-f", "--file"
            help = "Variant Consensus Fasta file"
            arg_type = String
            required = true
            "--variant-name"
            help = "Variant name"
            arg_type = String
            required = true
            "-w", "--window"
            help = "Slide window percent size to apply"
            arg_type = Float16
            required = true
            default = 0.01
            "-p", "--val-positions-file"
            help = "Positions file for validation"
            arg_type = String
            required = true
            "-o", "--output"
            help = "Output directory"
            arg_type = String
            required = true
        end

        parsed_args = parse_args(ARGS, setting)

        filePath::String = parsed_args["file"]
        windowSize::Float16 = parsed_args["window"]
        positionsBedFile::String = parsed_args["val-positions-file"]
        outputFile::String = parsed_args["output"]
        varname::String = parsed_args["variant-name"]

        println("Parsed args:")
        for (arg, val) in parsed_args
            println("  $arg  =>  $val")
        end

        positions::Vector{Int} = DataIO.readVectorFromFile(positionsBedFile, Int)
        validateEntropyWindow(positions, windowSize, varname, outputFile, filePath)
    end

    validate_region_main()

end

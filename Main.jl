
using Distributed, Pkg

include("DataIO.jl")
include("KmerUtils.jl")
include("EntropyUtil.jl")
include("TransformUtils.jl")

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

begin

    function _extractFreqWindow(
        numSeries::Array{Vector{T}},
        seqLen::Int
    )::Tuple{Int,Int} where {T<:Real}

        freqIndexes = Integer[]
        crossEspectrum = TransformUtils.elementWiseMult(numSeries, seqLen)
        N = MinMax(crossEspectrum)
        norm = N(crossEspectrum)

        @inbounds for ii in eachindex(norm)
            if norm[ii] >= 0.1
                push!(freqIndexes, ii)
            else
                crossEspectrum[ii] = zero(T)
            end
        end

        # normal = irfft(crossEspectrum, seqLen - 1)

        # plt = plot([norm, normal], layout=(2, 1))
        # savefig(plt, "filter2.png")

        # Filtro passa faixa, ou retornando apenas a frequencia como impulso
        if length(freqIndexes) > 0
            if (freqIndexes[1] == freqIndexes[length(freqIndexes)])
                return (1, freqIndexes[1])
            end
            return (freqIndexes[1], freqIndexes[length(freqIndexes)])
        else
            return (1, seqLen)
        end
    end


    function filterSignal!(
        points::Vector{Int},
        wndw::Tuple{Int,Int},
        wndwSize::Int,
        sequence::Vector{Float64},
        threshold::T
    ) where {T<:Real}

        fft = abs.(rfft(sequence))
        for ii in eachindex(fft)
            if (ii < wndw[1] || ii > wndw[2])
                fft[ii] = zero(T)
            end
        end
        normal = irfft(fft, length(sequence))
        N = MinMax(normal)
        norm = N(normal)

        # plt = plot([sequence, fft, norm], layout=(3, 1))
        # savefig(plt, "filter_current.png")

        for ii in eachindex(norm)
            if norm[ii] >= threshold
                push!(points, ii + wndwSize)
            end
        end

    end

    function extractRanges(norm::Vector{Int}, tolerance::Int)::Vector{Tuple{Int,Int}}
        ranges = Vector{Tuple{Int,Int}}()
        n = length(norm)
        if n == 0
            return ranges
        end
        start_idx = 1
        for i in 2:n
            # Check if the current value is not within the tolerance to the previous one
            if norm[i] > norm[i-1] + tolerance + 1
                # Add the current range (start_idx to i-1)
                push!(ranges, (norm[start_idx], norm[i-1]))
                # Start a new range from the current index
                start_idx = i
            end
        end
        push!(ranges, (norm[start_idx], norm[n]))
        return ranges
    end



    #WARNING: tem que considerar os finais das sequencias nos indices que não completam
    function _runMethodology!(
        filePath::String,
        slidWindw::Int,
        tolerance::Int,
        signalThreshold::Float64
    )

        seqStat = DataIO.getShortestLength(filePath)
        minSize = seqStat[1]
        initI = 1

        discPoints = Int[]
        println("Processing Regions")
        toCross = Array{Vector{Float64}}(undef, seqStat[2])
        @inbounds while (initI + slidWindw) <= minSize
            endI = (initI + slidWindw) - 1
            for (seqno, record) in enumerate(open(FASTAReader, filePath))
                toCross[seqno] = DataIO.sequence2NumericalSerie(sequence(record), initI, endI)
            end

            freqWindow = _extractFreqWindow(toCross, slidWindw)
            for serie in toCross
                filterSignal!(discPoints, freqWindow, initI - 1, serie, signalThreshold)
            end
            DataIO.progressBar!(endI, minSize)
            initI = endI + 1
        end
        DataIO.progressBar!(100, 100)

        ranges = extractRanges(sort(unique(discPoints)), tolerance)
        # @show ranges
        DataIO.writeFASTASingleChr!(filePath, "_output.fasta", ranges)
        println("")

    end

    function _runMethodology!(
        dirPath::String,
        regions::Vector{Tuple{String,Int,Int}},
        tolerance::Int,
        signalThreshold::Float64
    )
        files::Vector{String} = readdir(dirPath)

        extractRegionPoints = Dict{String,Vector{Tuple{Int,Int}}}()
        println("Applying RRM at regions")
        # Mudar no futuro para parelelizar so é possivel atualmente pq está linear
        toCross = Array{Vector{Float64}}(undef, length(files))
        for (i, (chr, initI, endI)) in enumerate(regions)
            DataIO.progressBar!(i, length(regions))
            discPoints = Int[]
            # Each file is a genome from a organism
            for (fileno, file) in enumerate(files)
                reader = open(FASTAReader, dirPath * "/" * file)
                for record in reader
                    if (identifier(record) == chr)
                        toCross[fileno] = DataIO.sequence2NumericalSerie(sequence(record), initI, endI)
                    end
                end
                close(reader)
            end

            freqWindow = _extractFreqWindow(toCross, endI - initI)
            for serie in toCross
                filterSignal!(discPoints, freqWindow, initI - 1, serie, signalThreshold)
            end
            ranges = extractRanges(sort(unique(discPoints)), tolerance)

            try
                extractRegionPoints[chr] = vcat(extractRegionPoints[chr], ranges)
            catch
                extractRegionPoints[chr] = ranges
            end
        end
        DataIO.progressBar!(1, 1)
        println("")
        DataIO.writeFASTAS!(dirPath, "_output.fasta", extractRegionPoints)

    end


    function main()
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
        println(parsed_args)

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
            # _runMethodology!(filePath, sldWindow, tolerance, threshold)
        else
            regions = DataIO.readRegionsFromBed(regionsbedfile)
            _runMethodology!(dirPath, regions, tolerance, threshold)
        end
    end

    function runEntropy()
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


    runEntropy()


end

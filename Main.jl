
using Distributed, Pkg
Pkg.instantiate()
Pkg.activate(".")

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
        testLbl::String,
        variant::String
    )
        dataset::String = "../datasets/consensus"

        # files::Vector{String} = readdir(dataset)
        # histograms = Vector{Vector{Float64}}()

        # para cada amostra da variante
        # for (i, fastaFile) in enumerate(files)
        @show variant

        sequences::Array{String} = []
        for record in open(FASTAReader, "$dataset/$variant")
            seq::String = sequence(String, record)
            push!(sequences, seq)
        end
        wndwStep::Int8 = 1
        plt = plot(title="Regions - $(wnwPercent*100)%")

        # Processo para encontrar valores de entropia por região do genoma
        for seqs in sequences
            slideWndw::Int = ceil(Int, length(seqs) * wnwPercent)
            y::Vector{Float64} = mountEntropyByWndw(slideWndw, wndwStep, seqs)
            N = MinMax(y)
            norm = N(y)
            # push!(histograms, N(y))

            # histograms[i] = N(y)
            x = range(1, length(norm))
            lim = [0, length(norm)]
            regions::Vector{Int} = histogramPosWndw(positions, slideWndw, length(norm))

            plot!(twinx(), x, regions,
                label="Frequency",
                seriestype=:bar,
                linecolor=nothing,
                xlims=lim)
            plot!(x, norm, label="Entropy-value", xlims=lim)
            # plot!(x, findPosWndw(positions, slideWndw, norm), label="Pos-existence", xlims=lim)
        end

        png(plt, testLbl)
        # end
    end

    # GAMMA VARIANT ANNOTATION ASSERT
    gammaPositions_90 = Vector{Int}([241,
        733,
        2749,
        3037,
        5648,
        6319,
        6613,
        11287,
        12778,
        13860,
        14408,
        17259,
        21614,
        21621,
        21638,
        21974,
        22132,
        22812,
        23012,
        23063,
        23403,
        23525,
        24642,
        25088,
        26149,
        28167,
        28512,
        28877,
        28878,
        28881,
        28882,
        28883])
    gammaPositions_99 = Vector{Int}([3037, 5648, 14408, 23403, 23525, 25088, 26149, 28512])

    # ALPHA VARIAN ASSERTION
    alphaPositions_90 = Vector{Int}([241,
        2470,
        2832,
        3037,
        5386,
        8393,
        10029,
        10449,
        11537,
        13195,
        14408,
        15240,
        18163,
        23202,
        23403,
        23525,
        23599,
        23604,
        23854,
        23948,
        24130,
        24424,
        24469,
        24503,
        25584,
        26270,
        26577,
        26709,
        27259,
        28881,
        28882,
        28883])
    alphaPositions_99 = Vector{Int}([10029, 23403, 3037])

    # OMICRON VARIANT ANNOTATION
    omicronPositions_90 = Vector{Int}([241,
        670,
        2790,
        2832,
        3037,
        4184,
        4321,
        5386,
        6512,
        8393,
        9344,
        9424,
        9534,
        9866,
        10029,
        10198,
        10447,
        10449,
        11282,
        11287,
        11537,
        12880,
        13195,
        14408,
        15240,
        15714,
        17410,
        18163,
        19955,
        20055,
        21618,
        21762,
        21846,
        21987,
        22200,
        22578,
        22674,
        22679,
        22686,
        22688,
        23063,
        23075,
        23202,
        23403,
        23525,
        23599,
        23604,
        23854,
        23948,
        24130,
        24424,
        24469,
        24503,
        25000,
        25584,
        26060,
        26270,
        26577,
        26709,
        26858,
        27259,
        27382,
        27383,
        27384,
        27807,
        28271,
        28311,
        28361,
        28881,
        28882,
        28883,
        29510])
    omicronPositions_99 = Vector{Int}([
        2790,
        3037,
        5386,
        8393,
        10029,
        12880,
        13195,
        14408,
        15714,
        17410,
        23403,
        23525,
        24424,
        24469,
        26060])

    windows = Vector{Float16}([0.1, 0.15, 0.2, 0.25, 0.3])

    for w in windows
        validateEntropyWindow(gammaPositions_90, w, "Gamma-histogram90-wndn=$w", "Gamma/Gamma_reference.fasta")
        validateEntropyWindow(gammaPositions_99, w, "Gamma-histogram99-wndn=$w", "Gamma/Gamma_reference.fasta")

        validateEntropyWindow(alphaPositions_90, w, "Alpha-histogram90-wndn=$w", "Alpha/Alpha_reference.fasta")
        validateEntropyWindow(alphaPositions_99, w, "Alpha-histogram99-wndn=$w", "Alpha/Alpha_reference.fasta")

        validateEntropyWindow(omicronPositions_90, w, "Omicron-histogram90-wndn=$w", "Omicron/Omicron_reference.fasta")
        validateEntropyWindow(omicronPositions_99, w, "Omicron-histogram99-wndn=$w", "Omicron/Omicron_reference.fasta")
    end


end

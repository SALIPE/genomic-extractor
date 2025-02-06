module RRM
include("DataIO.jl")
include("TransformUtils.jl")

using .DataIO,
    .TransformUtils,
    Normalization,
    AbstractFFTs,
    Plots

export RRM

function getFourierCoefficient(
    cuttedSequences::Array{String},
    seqLen::Int64
)

    toCross = Array{Vector{Float64}}(undef, length(cuttedSequences))

    for (seqno, sequence) in enumerate(cuttedSequences)
        toCross[seqno] = DataIO.sequence2NumericalSerie(sequence)
    end

    crossEspectrum = TransformUtils.elementWiseMult(toCross, seqLen)

    return crossEspectrum

end

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

# ------------ RRM ---------------
# Funções que rodam slide window para executar o RRM spectrum 


function RRMEntropySignal(
    consensusSignals::Vector{Tuple{String,Vector{Float64}}}
)
    signals = Vector{Tuple{String,Vector{Float64}}}(undef, length(consensusSignals))

    # plt = plot(title="Entropy Signals Norms.", dpi=300)
    for (i, (class, signal)) in enumerate(consensusSignals)
        N = MinMax(signal)
        norm = N(signal)
        # plot!(norm, label=class)
        signals[i] = (class, norm)
    end
    # png(plt, "norm_signals")

    min_length = minimum(map(x -> length(x[2]), signals))

    freqWindow = _extractFreqWindow(map(x -> x[2], signals), min_length)

    for (i, (class, sequence)) in enumerate(signals)
        fft = abs.(rfft(sequence))
        for ii in eachindex(fft)
            if (ii < freqWindow[1] || ii > freqWindow[2])
                fft[ii] = zero(fft[ii])
            end
        end
        normal = irfft(fft, length(sequence))
        signals[i] = (class, normal)
    end

    return signals

end

#WARNING: tem que considerar os finais das sequencias nos indices que não completam
function runRRMMethodology!(
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
    while (initI + slidWindw) <= minSize
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

function runRRMMethodology!(
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
end
# -------------------------------
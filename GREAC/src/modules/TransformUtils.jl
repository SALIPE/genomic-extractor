module TransformUtils

include("DataIO.jl")

using .DataIO, AbstractFFTs, Normalization

export TransformUtils


function elementWiseMult(
    series::Array{Vector{T}}
)::Vector{T} where {T<:Real}

    crossEspectrum::Vector{T} = abs.(rfft(series[1])[2:end])
    for s in series[2:end]
        dft = abs.(rfft(s)[2:end])
        crossEspectrum = crossEspectrum .* dft
    end
    return crossEspectrum
end

function getFourierCoefficient(
    cuttedSequences::Array{String}
)::Vector{Float64}

    toCross = Array{Vector{Float64}}(undef, length(cuttedSequences))

    @inbounds for (i, sequence) in enumerate(cuttedSequences)
        toCross[i] = DataIO.sequence2AminNumSerie(sequence)
    end

    crossEspectrum = TransformUtils.elementWiseMult(toCross)
    N = MinMax(crossEspectrum)
    norm = N(crossEspectrum)

    return norm

end

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

function _extractFreqWindow(
    numSeries::Array{Vector{T}},
    seqLen::Int
)::Tuple{Int,Int} where {T<:Real}

    freqIndexes = Integer[]
    crossEspectrum = elementWiseMult(numSeries)
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

end
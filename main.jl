include("dataIO.jl")
include("transformUtils.jl")

using AbstractFFTs, FASTX, Plots, LoopVectorization, Normalization
using .DataIO, .TransformUtils


function _extractFreqWindow(
    numSeries::Array{Vector{T}},
    seqLen::Int
)::Vector{Int} where {T<:Real}

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

    normal = irfft(crossEspectrum, seqLen - 1)

    plt = plot([norm, normal], layout=(2, 1))
    savefig(plt, "filter2.png")

    if (freqIndexes[1] == freqIndexes[length(freqIndexes)])
        return [1, freqIndexes[1]]
    end
    return [freqIndexes[1], freqIndexes[length(freqIndexes)]]

end

let
    filePath::String = "/home/salipe/Desktop/GitHub/datasets/gen_dron_car.fasta"
    seqReader = open(FASTAReader, filePath)


    # minLength = DataIO.getShortestLength(filePath)
    powder = 512

    initI = 10000
    endI = initI + powder

    toCross = Array{Vector{Float64}}(undef, 3)
    for (seqno, record) in enumerate(seqReader)
        toCross[seqno] = DataIO.sequence2NumericalSerie(sequence(record), initI, endI)
    end

    freqWindow = _extractFreqWindow(toCross, powder)
    @show freqWindow

    # numSeries = DataIO.sequence2NumericalSerie(seq1)
    # dft::Vector{ComplexF64} = rfft(numSeries)

    # filtered = zeros(ComplexF64, length(dft))
    # ii = firstindex(dft)
    # @inbounds while ii <= 25
    #     filtered[ii] = dft[ii]
    #     ii += 1
    # end
    # normal = irfft(filtered, length(numSeries))
    # N = MinMax(normal)
    # cross = N(normal)
    # plt = plot([numSeries, abs.(dft[2:length(dft)]), abs.(filtered[2:length(dft)]), cross], layout=(4, 1))
    # savefig(plt, "filter.png")

    # idx_commom = Vector{Int}()
    # i = firstindex(numSeries)
    # @inbounds while i <= lastindex(numSeries)
    #     rounded = round(normal[i], digits=4)
    #     if numSeries[i] == rounded
    #         push!(idx_commom, i)
    #     end
    #     i += 1
    # end
    # @show idx_commom

end


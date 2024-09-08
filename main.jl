include("dataIO.jl")
include("transformUtils.jl")

using AbstractFFTs, FASTX, Plots, LoopVectorization, Normalization
using .DataIO, .TransformUtils



function _separateSequences!(
    series::Array{Vector{Float64}}
)::Vector{Float64}

    minLength::Int = length(series[argmin(series)])

    toCross::Array{Vector{Float64}} = []
    toReprocess::Array{Vector{Float64}} = []
    for numSeries in series
        seqLen::Int = length(seq)
        if seqLen == minLength
            push!(toCross, numSeries)
        else
            push!(crossParts, numSeries[1:minLength])
            push!(toReprocess, numSeries[minLength:length(numSeries)])
        end
    end

    if length(toReprocess) > 1
        return _separateSequences!(toReprocess)
    end

    # cross = TransformUtils.elementWiseMult(toCross)

end




let
    # filePath::String = "/home/salipe/Desktop/GitHub/datasets/gen_dron_car.fasta"
    # sequences = open(FASTAReader, filePath)

    # minLength = DataIO.getShortestLength(filePath)
    powder = 256

    seq1::String = "GTATCCACAAAGTTAATTATCGCGAACGAGTCCTCCTTTCAACCGTGCCACGCGTGCCGCACTCGATGCGAGGAAGAAATTTCCGTTTTCAAAAGCCGTCTGTCTCTCCATTATCCGCTACCATCCTCTCGTACCACCATCTCTTTCTCGTGCCCCAAATCACCCC"
    numSeries = DataIO.sequence2NumericalSerie(seq1)
    dft::Vector{ComplexF64} = rfft(numSeries)

    filtered = zeros(ComplexF64, length(dft))
    ii = firstindex(dft)
    @inbounds while ii <= 25
        filtered[ii] = dft[ii]
        ii += 1
    end
    normal = irfft(filtered, length(numSeries))
    N = MinMax(normal)
    cross = N(normal)
    plt = plot([numSeries, abs.(dft[2:length(dft)]), abs.(filtered[2:length(dft)]), cross], layout=(4, 1))
    savefig(plt, "filter.png")

    idx_commom = Vector{Int}()
    i = firstindex(numSeries)
    @inbounds while i <= lastindex(numSeries)
        rounded = round(normal[i], digits=4)
        if numSeries[i] == rounded
            push!(idx_commom, i)
        end
        i += 1
    end
    @show idx_commom

    # toCross::Array{Vector{Float64}} = []
    # toReprocess::Array{Vector{Float64}} = []
    # for seq in sequences
    #     seqLen::Int = seqsize(seq)
    #     numSeries = DataIO.sequence2NumericalSerie(sequence(seq))
    #     push!(toCross, numSeries[25000000:25000000+powder])
    #     # if seqLen == minLength
    #     #     push!(toCross, numSeries)
    #     # else
    #     #     push!(toCross, numSeries[1:minLength])
    #     #     push!(toReprocess, numSeries[minLength:length(numSeries)])
    #     # end
    # end
    # cross = TransformUtils.elementWiseMult(toCross, powder)
    # N = MinMax(cross)
    # cross = N(cross)
    # dftfreq = rfftfreq(minLength)

    # plt = plot(cross, xlabel="Frequency (Hz)", ylabel="Magnitude", title="FFT of the Signal")


    # if length(toReprocess) > 1
    #     return _separateSequences!(toReprocess)
    # end


end


# seqLen = length(convsig)
# less::Int8 = seqLen % 2 == 0 ? 1 : 2

# dftfreq = rfftfreq(2000)
# plt = plot(dftfreq, abs.(dftv[1000:2000]), xlabel="Frequency (Hz)", ylabel="Magnitude", title="FFT of the Signal")
# savefig(plt, "myplotconv2.png")


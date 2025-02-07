module Classification

include("DataIO.jl")

using .DataIO,
    FLoops,
    Serialization

export Classification

function classifyInput(
    inputSequence::AbstractString,
    outputfilename::Union{Nothing,AbstractString}=nothing
)

    @show outputfilename
    modelCachedFile = "$(pwd())/.project_cache/trained_model.dat"

    model::Union{Nothing,Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}} = DataIO.load_cache(modelCachedFile)

    if !isnothing(model)
        @info "Using model from cached data from $modelCachedFile"
        classifyInput(inputSequence, model, outputfilename)
    else
        error("Model not found in cached files!")
    end


end

function classifyInput(
    inputSequence::Base.CodeUnits,
    #Dict{VariantName, (Marked, Fourier Coefficients, Kmers)}
    model::Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}},
    outputfilename::Union{Nothing,String}
)

    report = Dict{String,Vector{Tuple{Tuple{UInt16,UInt16},UInt16}}}()
    for (key, (marked, _, kmers)) in model
        inputlen = minimum(length, [inputSequence, marked])
        report[key] = Vector{Tuple{Tuple{UInt16,UInt16},UInt16}}()
        limitedMark::BitArray = marked[1:inputlen]
        start = 0
        current = false

        for (i, bit) in enumerate(limitedMark)
            if bit && !current
                start = i
                current = true
            elseif !bit && current
                current = false
                count = @views countPatterns(inputSequence[start:i-1], kmers)

                push!(report[key], ((start, i - 1), count))
            end
        end
        if current
            count = @views countPatterns(inputSequence[start:inputlen], kmers)
            push!(report[key], ((start, inputlen), count))
        end
    end

    reportFilename::String = isnothing(outputfilename) ? "report" : outputfilename

    open("$(pwd())/$(reportFilename).txt", "w") do file

        for (var, regions) in report
            write(file, "\n\n########### $(uppercase(var)) ############")
            write(file, "\nTotal Exclusive Kmers: $(length(model[var][3]))")
            write(file, "\n####################################\n")
            for ((initPos, endPos), count) in regions
                write(file, "\nWindow Position( $initPos - $endPos ): \nKmer Count: $count")
            end
        end
    end

end

function countPatterns(
    seqWindow::SubArray,
    kmers::Vector{String})::UInt8

    patterns = [Base.Fix1(occursinKmerBit, codeunits(kmer)) for kmer in kmers]
    count::UInt8 = 0

    @floop for pattern in patterns
        if pattern(seqWindow)
            @reduce count += 1
        end
    end
    return count
end

function occursinKmerBit(
    kmer::Base.CodeUnits,
    windowBuffer::SubArray
)::Bool
    wlen = length(windowBuffer)
    klen = length(kmer)

    @inbounds for i in 1:(wlen-klen+1)
        match = true
        for j in 1:klen
            (windowBuffer[i+j-1] ≠ kmer[j]) && (match = false; break)
        end
        match && return true
    end
    return false

end

end

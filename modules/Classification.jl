module Classification

include("DataIO.jl")

using .DataIO,
    FLoops,
    Serialization

export Classification

function classifyInput(
    inputFile::AbstractString,
    exclusiveKmers::Dict{String,Vector{String}},
)
    modelCachedFile = "$(pwd())/.cache/trained_model.dat"

    model::Union{Nothing,} = DataIO.load_cache(modelCachedFile)

    if !isnothing(cache)
        @info "Using model from cached data from $cache_path"
        return classifyInput(inputFile, exclusiveKmers, model)
    else
        error("Model not found in cached files!")
    end


end

function classifyInput(
    inputSequence::AbstractString,
    exclusiveKmers::Dict{String,Vector{String}},
    #Dict{VariantName, (Marked, Fourier Coefficients)}
    model::Dict{String,Tuple{BitArray,Vector{Vector{Float64}}}}
)

    report = Dict{String,Vector{Tuple{Tuple{Uint16,Uint16},UInt16}}}()
    for (key, (marked, coefs)) in model
        inputlen = minimum(length, [inputSequence, marked])

        limitedMark::BitArray = marked[1:inputlen]
        start = 0
        current = false

        for (i, bit) in enumerate(limitedMark)
            if bit && !current
                start = i
                current = true
            elseif !bit && current
                current = false

                count = countPatterns(codeunits(inputSequence[start:i-1]), exclusiveKmers[key])

                push!(report[key], ((start, i - 1), count))
            end
        end
        if current
            count = countPatterns(codeunits(inputSequence[start:inputlen]), exclusiveKmers[key])
            push!(report[key], ((start, inputlen), count))
        end
    end


    open("$(pwd())/report.txt", "w") do file

        for (var, regions) in report
            write(file, "\nVariant: $var")

            for (i, (_, count)) in enumerate(regions)
                write(file, "\n$i Region: $count")
            end
        end
    end



end

function countPatterns(
    seqWindow::Base.CodeUnits,
    kmers::Vector{String})::UInt8

    patterns = [Base.Fix1(occursinKmerBit, (codeunits(kmer), length(seqWindow))) for kmer in kmers]

    count::UInt8 = 0

    @floop for pattern in patterns
        if pattern(seqWindow[initPos:endPos])
            @reduce count += 1
        end
    end
    return count
end

function occursinKmerBit(
    args,
    window_buffer
)::Bool
    kmer = args[1]
    wndwSize = args[2]
    klen = length(kmer)
    # Slide through window to find matches
    for pos in 1:(wndwSize-klen+1)
        match = true
        for i in 1:klen
            if window_buffer[pos+i-1] ≠ kmer[i]
                match = false
                break
            end
        end
        if match
            return true
        end
    end
    return false
end

end
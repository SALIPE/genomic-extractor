module Classification

include("DataIO.jl")

using .DataIO,
    FLoops,
    Serialization

export Classification

function classifyInput(
    inputFile::AbstractString,
)
    modelCachedFile = "$(pwd())/.cache/trained_model.dat"

    model::Union{Nothing,Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}} = DataIO.load_cache(modelCachedFile)

    if !isnothing(model)
        @info "Using model from cached data from $cache_path"
        return classifyInput(inputFile, model)
    else
        error("Model not found in cached files!")
    end


end

function classifyInput(
    inputSequence::AbstractString,
    #Dict{VariantName, (Marked, Fourier Coefficients, Kmers)}
    model::Dict{String,Tuple{BitArray,Vector{Vector{Float64}},Vector{String}}}
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

                count = countPatterns(codeunits(inputSequence[start:i-1]), kmers)

                push!(report[key], ((start, i - 1), count))
            end
        end
        if current
            count = countPatterns(codeunits(inputSequence[start:inputlen]), kmers)
            push!(report[key], ((start, inputlen), count))
        end
    end

    open("$(pwd())/report.txt", "w") do file

        for (var, regions) in report
            write(file, "\nVariant: $var - Exclusive Kmers")

            for ((initPos, endPos), count) in regions
                write(file, "\nWindow Position( $initPos - $endPos ): \nKmer Count: $count")
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
        if pattern(seqWindow)
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

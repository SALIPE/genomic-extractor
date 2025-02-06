module Classification

include("DataIO.jl")

using .DataIO, Serialization

export Classification

function classifyInput(
    inputFile::AbstractString
)
    modelCachedFile = "$(pwd())/.cache/trained_model.dat"

    model::Union{Nothing,} = DataIO.load_cache(modelCachedFile)

    if !isnothing(cache)
        @info "Using model from cached data from $cache_path"
        return classifyInput(inputFile, model)
    else
        error("Model not found in cached files!")
    end


end

function classifyInput(
    inputSequence::AbstractString,
    #Dict{VariantName, (Marked, Fourier Coefficients)}
    model::Dict{String,Tuple{BitArray,Vector{Vector{Float64}}}}
)

    report = Dict{String,Tuple{UInt16,UInt16}}()
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
            end
        end
    end

end

end
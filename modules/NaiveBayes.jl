module NaiveBayes

include("Model.jl")

using FLoops, .Model
export NaiveBayes

struct MultiClassNaiveBayes
    classes::Vector{String}
    priors::Dict{String,Float64}
    class_string_probs::Dict{String,Vector{Int}}
    wnw_size::Int
    max_seq_windows::Int
    kmerset::Set{String}
end

function fitMulticlassNB(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int,
    max_seq_windows::Int
)::MultiClassNaiveBayes

    priors = Dict{String,Float64}()
    class_string_probs = Dict{String,Vector{Int}}()

    total_samples = sum(x -> x[2], meta_data)

    for (class, seq_total) in meta_data

        println("Calculating $class probabilities")

        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (wnw_size, max_seq_windows, byte_seqs[class]))

        @floop for kmer in collect(kmerset)
            kmer_seq_histogram = get_class_appearences(kmer)

            @reduce(
                kmer_distribution = zeros(UInt64, max_seq_windows) .+ kmer_seq_histogram
            )
        end

        class_string_probs[class] = kmer_distribution

        # Process Overall Frequence
        # class_string_probs[class] = kmer_distribution ./ (length(kmerset) * length(byte_seqs[class]))
        priors[class] = seq_total / total_samples
    end

    return MultiClassNaiveBayes(
        [class for (class, _) in meta_data],
        priors,
        class_string_probs,
        wnw_size,
        max_seq_windows,
        kmerset)
end



function def_kmer_classes_probs(
    seq_data::Tuple{Int,Int,Vector{Base.CodeUnits}},
    kmer::String)::Vector{UInt64}

    wnw_size, max_seq_windows, sequences = seq_data

    fn_occursin = Base.Fix1(Model.occursinKmerBit, codeunits(kmer))

    @floop for seq in sequences
        seq_windows = length(seq) - wnw_size + 1
        local_seq_histogram = zeros(UInt64, max_seq_windows)
        for initPos in 1:seq_windows
            endPos = initPos + wnw_size - 1
            wndw_buffer = @view seq[initPos:endPos]

            if fn_occursin(wndw_buffer)
                local_seq_histogram[initPos] += 1
            end
        end

        @reduce(
            seq_histogram = zeros(UInt64, max_seq_windows) .+ local_seq_histogram
        )

    end

    return seq_histogram
end

function def_kmer_presence(
    seq_data::Tuple{Int,Int,Base.CodeUnits},
    kmer::String
)::Tuple{String,BitArray}

    wnw_size, max_seq_windows, seq = seq_data

    fn_occursin = Base.Fix1(Model.occursinKmerBit, codeunits(kmer))
    seq_presence = falses(max_seq_windows)

    seq_windows = length(seq) - wnw_size + 1

    for initPos in 1:seq_windows
        endPos = initPos + wnw_size - 1
        wndw_buffer = @view seq[initPos:endPos]

        if fn_occursin(wndw_buffer)
            seq_presence[initPos] = 1
        end
    end

    return (kmer, seq_presence)

end


function predict(
    model::MultiClassNaiveBayes,
    X::Vector{Float64}
)::Tuple{String,Dict{String,Float64}}
    log_probs = Dict{String,Float64}([(class, zero(Float64)) for class in model.classes])


    for c in model.classes
        log_prob = log(model.priors[c])

        # Get the class's precomputed probabilities for this string
        class_probs = model.class_string_probs[c]

        # Compute the dot product between input probabilities and class probabilities
        likelihood = sum(X .* class_probs)
        # Avoid log(0) by adding a small epsilon (e.g., 1e-9)
        likelihood = max(likelihood, 1e-9)

        log_prob += log(likelihood)

        log_probs[c] = log_prob
    end

    # Return the class with the highest log probability
    return argmax(log_probs), log_probs
    # return log_probs
end




end
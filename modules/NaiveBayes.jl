module NaiveBayes

include("Model.jl")

using FLoops, .Model
export NaiveBayes

struct MultiClassNaiveBayes
    classes::Vector{String}
    priors::Dict{String,Float64}
    class_string_probs::Dict{String,Dict{String,Vector{Float64}}}
end

function fitMulticlassNB(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int
)::MultiClassNaiveBayes


    seq_windows = maximum(x -> length(x[2]), byte_seqs) - wnw_size + 1

    @info seq_windows
    priors = Dict{String,Float64}()
    class_string_probs = Dict{String,Dict{String,Vector{Float64}}}()

    total_samples = sum(x -> x[2], meta_data)

    for (class, seq_total) in meta_data
        get_class_appearences = Base.Fix1(def_kmer_classes_probs, byte_seqs[class])

        create_kmer_distribution = [Base.Fix1(get_class_appearences, kmer) for kmer in collect(kmerset)]

        kmer_distribution = Dict{String,Vector{Float64}}([(kmer, Vector{Float64}(undef, seq_windows)) for kmer in kmerset])
        @floop for get_appearences in create_kmer_distribution
            kmer, seq_histogram = get_appearences()
            kmer_distribution[kmer] = seq_histogram ./ seq_total
        end
        class_string_probs[class] = seq_total / total_samples
    end

    return MultiClassNaiveBayes(keys(meta_data), priors, class_string_probs)
end

function def_kmer_classes_probs(
    seq_data::Tuple{Int,Vector{Base.CodeUnits}},
    kmer::String)::Tuple{String,Vector{UInt64}}

    seq_windows, sequences = seq_data

    fn_occursin = Base.Fix1(Model.occursinKmerBit, codeunits(kmer))
    seq_histogram = zeros(UInt64, seq_windows)

    for seq in sequences
        seq_windows = length(seq) - data.wnw_size + 1

        for initPos in 1:seq_windows
            endPos = initPos + data.wnw_size - 1
            wndw_buffer = @view seq[initPos:endPos]

            if fn_occursin(wndw_buffer)
                seq_histogram[initPos] += 1
            end
        end

    end

    return (kmer, seq_histogram)
end


function predict(
    model::MultiClassNaiveBayes,
    X::Dict{String,Vector{Float64}}  # Input: Probabilities for each string
)
    log_probs = Dict{Any,Float64}()

    for c in model.classes
        log_prob = log(model.priors[c])  # Log prior

        for (string, input_probs) in X
            # Get the class's precomputed probabilities for this string
            class_probs = model.class_string_probs[c][string]

            # Compute the dot product between input probabilities and class probabilities
            likelihood = sum(input_probs .* class_probs)

            # Avoid log(0) by adding a small epsilon (e.g., 1e-9)
            likelihood = max(likelihood, 1e-9)

            log_prob += log(likelihood)
        end

        log_probs[c] = log_prob
    end

    # Return the class with the highest log probability
    return argmax(log_probs)
end






end
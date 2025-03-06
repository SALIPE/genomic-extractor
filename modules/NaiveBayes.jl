module NaiveBayes

include("Model.jl")

using FLoops, .Model
export NaiveBayes

struct MultiClassNaiveBayes
    classes::Vector{String}
    priors::Dict{String,Float64}
    class_string_probs::Dict{String,Dict{String,Vector{Float64}}}
    wnw_size::Int
    max_seq_windows::Int
end

function fitMulticlassNB(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int,
    max_seq_windows::Int
)::MultiClassNaiveBayes


    priors = Dict{String,Float64}()
    class_string_probs = Dict{String,Dict{String,Vector{Float64}}}()

    total_samples = sum(x -> x[2], meta_data)

    for (class, seq_total) in meta_data

        println("Calculating $class probabilities")

        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (wnw_size, max_seq_windows, byte_seqs[class]))

        kmer_distribution = Dict{String,Vector{Float64}}([(kmer, Vector{Float64}(undef, max_seq_windows)) for kmer in kmerset])

        @floop for kmer in collect(kmerset)
            kmer, seq_histogram = get_class_appearences(kmer)
            kmer_distribution[kmer] = seq_histogram ./ seq_total
        end

        class_string_probs[class] = kmer_distribution
        priors[class] = seq_total / total_samples
    end

    return MultiClassNaiveBayes([class for (class, _) in meta_data], priors, class_string_probs, wnw_size, max_seq_windows)
end

function def_kmer_classes_probs(
    seq_data::Tuple{Int,Int,Vector{Base.CodeUnits}},
    kmer::String)::Tuple{String,Vector{UInt64}}

    wnw_size, max_seq_windows, sequences = seq_data

    fn_occursin = Base.Fix1(Model.occursinKmerBit, codeunits(kmer))
    seq_histogram = zeros(UInt64, max_seq_windows)

    @floop for seq in sequences
        seq_windows = length(seq) - wnw_size + 1

        for initPos in 1:seq_windows
            endPos = initPos + wnw_size - 1
            wndw_buffer = @view seq[initPos:endPos]

            if fn_occursin(wndw_buffer)
                seq_histogram[initPos] += 1
            end
        end

    end

    return (kmer, seq_histogram)
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
    X::Dict{String,BitArray}
)::Dict{String,Float64}
    log_probs = Dict{String,Float64}([(class, zero(Float64)) for class in model.classes])

    likelihood_fns = [Base.Fix1(getLikelihood, (string, input_probs)) for (string, input_probs) in X]

    for c in model.classes
        log_prob = log(model.priors[c])

        @floop for get_likelihood in likelihood_fns
            likelihood = get_likelihood(model.class_string_probs[c])
            @reduce log_prob += log(likelihood)
        end

        log_probs[c] = log_prob
    end

    # Return the class with the highest log probability
    # return argmax(log_probs), log_probs
    return log_probs
end

function getLikelihood(
    input::Tuple{String,BitArray},
    class_string_probs
)
    (string, input_probs) = input
    # Get the class's precomputed probabilities for this string
    class_probs = class_string_probs[string]

    # Compute the dot product between input probabilities and class probabilities
    likelihood = sum(input_probs .* class_probs)
    # Avoid log(0) by adding a small epsilon (e.g., 1e-9)
    likelihood = max(likelihood, 1e-9)

    return likelihood
end


end
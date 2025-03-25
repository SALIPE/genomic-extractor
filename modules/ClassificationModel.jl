module ClassificationModel

include("RegionExtraction.jl")

using FLoops, .RegionExtraction, LinearAlgebra, Statistics
export ClassificationModel

struct MultiClassModel
    classes::Vector{String}
    priors::Dict{String,Float64}
    class_string_probs::Dict{String,Vector{Float64}}
    wnw_size::Int
    max_seq_windows::Int
    kmerset::Set{String}
    regions::Vector{Tuple{Int,Int}}
end


function fitMulticlass(
    kmerset::Set{String},
    kmers_dist::Dict{String,Int},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int,
    max_seq_windows::Int,
    regions::Vector{Tuple{Int,Int}}
)::MultiClassModel

    priors = Dict{String,Float64}()
    class_string_probs = Dict{String,Vector{Float64}}()

    total_samples = sum(x -> x[2], meta_data)
    regions_len = length(regions)
    for (class, _) in meta_data

        println("Calculating $class probabilities")

        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (regions, byte_seqs[class]))

        @floop for kmer in collect(kmerset)
            kmer_seq_histogram = get_class_appearences(kmer)

            @reduce(
                kmer_distribution = zeros(UInt64, regions_len) .+ kmer_seq_histogram
            )
        end

        # Process Overall Frequence
        class_string_probs[class] = kmer_distribution ./ (length(kmerset) * length(byte_seqs[class]))
        priors[class] = meta_data[class] / total_samples
    end

    return MultiClassModel(
        [class for (class, _) in meta_data],
        priors,
        class_string_probs,
        wnw_size,
        max_seq_windows,
        kmerset,
        regions)
end



function def_kmer_classes_probs(
    seq_data::Tuple{Vector{Tuple{Int,Int}},Vector{Base.CodeUnits}},
    kmer::String)::Vector{UInt64}

    regions, sequences = seq_data

    regions_len = length(regions)
    fn_occursin = Base.Fix1(RegionExtraction.occursinKmerBit, codeunits(kmer))

    @floop for seq in sequences
        local_seq_histogram = zeros(UInt64, regions_len)
        seq_len = length(seq)

        for i in eachindex(regions)
            init_pos, end_pos = regions[i]

            if (end_pos > seq_len)
                end_pos = seq_len
            end

            wndw_buffer = @view seq[init_pos:end_pos]

            if fn_occursin(wndw_buffer)
                local_seq_histogram[i] += 1
            end
        end

        @reduce(
            seq_histogram = zeros(UInt64, regions_len) .+ local_seq_histogram
        )

    end

    return seq_histogram
end

function def_kmer_presence(
    seq_data::Tuple{Int,Int,Base.CodeUnits},
    kmer::String
)::Tuple{String,BitArray}

    wnw_size, max_seq_windows, seq = seq_data

    fn_occursin = Base.Fix1(RegionExtraction.occursinKmerBit, codeunits(kmer))
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



function predict_raw(
    parameters::Tuple{MultiClassModel,Union{Nothing,String}},
    X::Vector{Float64},)::Tuple{String,Dict{String,Float64}}

    model, metric = parameters

    if (isnothing(metric))
        metric = "manhattan"
    end

    dists = Dict{String,Float64}([(class, zero(Float64)) for class in model.classes])
    epsilon = 1e-6

    for c in model.classes

        # Get the class's precomputed conditional frequencies
        class_freqs = model.class_string_probs[c]

        if metric == "manhattan"

            # Manhattan distance
            dists[c] = sum(abs.(X - class_freqs))

        elseif metric == "euclidian"
            # Euclidian distance
            dists[c] = sqrt(sum((X - class_freqs) .^ 2))

        elseif metric == "mahalanobis"

            train_data = hcat([model.class_string_probs[c] for c in model.classes]...)

            covariance = cov(train_data; dims=2)
            inv_covariance = inv(covariance + epsilon * I(size(covariance, 1)))

            # Mahalanobis distance requires inverse covariance matrix
            if inv_covariance === nothing
                error("Mahalanobis metric requires inverse covariance matrix")
            end

            delta = X - class_freqs
            dists[c] = sqrt(delta' * inv_covariance * delta)

        elseif metric == "chisquared"
            # Chi-squared distance
            dists[c] = sum((X - class_freqs) .^ 2 ./ (class_freqs .+ 1e-9))

        elseif metric == "kld"

            #  Kullback-Leibler (KL) divergence
            Q_norm = X ./ sum(X)
            P_norm = class_freqs ./ sum(class_freqs)

            # Smooth to avoid zeros
            P_smoothed = P_norm .+ epsilon
            P_smoothed = P_smoothed ./ sum(P_smoothed)

            # Compute KL(Q || P_smoothed)
            dists[c] = sum(q * (log(q) - log(p)) for (q, p) in zip(Q_norm, P_smoothed) if q > 0)
        else
            error("Unsupported metric: $metric")
        end

    end

    return argmin(dists), dists
end



end
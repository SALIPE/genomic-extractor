module ClassificationModel

include("RegionExtraction.jl")

using FLoops, .RegionExtraction, LinearAlgebra, Statistics, StatsBase
export ClassificationModel

struct MultiClassModel
    classes::Vector{String}
    priors::Dict{String,Float64}
    class_string_probs::Dict{String,Vector{Float64}}
    variant_stats::Dict{String,Dict{Symbol,Any}}
    wnw_size::Int
    kmerset::Set{String}
    regions::Vector{Tuple{Int,Int}}
end


function fitMulticlass(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int,
    regions::Vector{Tuple{Int,Int}}
)::MultiClassModel

    priors = Dict{String,Float64}()
    class_string_probs = Dict{String,Vector{Float64}}()
    variant_stats = Dict{String,Dict{Symbol,Any}}()


    total_samples = sum(x -> x[2], meta_data)
    regions_len = length(regions)
    for (class, _) in meta_data

        class_seqs::Vector{Base.CodeUnits} = byte_seqs[class]
        println("Calculating $class probabilities")
        get_class_appearences = Base.Fix1(def_kmer_classes_probs, (regions, class_seqs))

        @floop for kmer in collect(kmerset)
            kmer_seq_histogram = get_class_appearences(kmer)

            @reduce(
                kmer_distribution = zeros(UInt64, regions_len) .+ kmer_seq_histogram
            )
        end

        # Process F_wr
        class_freq = kmer_distribution ./ (length(kmerset) * length(class_seqs))
        class_string_probs[class] = class_freq
        priors[class] = meta_data[class] / total_samples

        in_group = Vector{Float64}(undef, length(class_seqs))
        @floop for s in eachindex(class_seqs)
            seq::Base.CodeUnits = class_seqs[s]

            seq_distribution = sequence_kmer_distribution(regions, seq, collect(kmerset)) ./ length(kmerset)

            in_group[s] = sum(abs.(seq_distribution - class_freq))

        end

        variant_stats[class] = Dict(
            :mu => mean(in_group),
            :sigma => std(in_group),
            :percentiles => quantile(in_group, [0.05, 0.95])
        )
    end

    return MultiClassModel(
        [class for (class, _) in meta_data],
        priors,
        class_string_probs,
        variant_stats,
        wnw_size,
        kmerset,
        regions)
end

function trapezoidal_membership(
    stats::Dict{Symbol,Any},
    d::Float64)

    p5, p95 = stats[:percentiles]

    if d <= p5
        0.0
    elseif p5 < d <= p95
        1.0
    else
        max(0.0, 1 - (d - p95) / (p95 - p5))
    end
end

function gaussian_membership(
    stats::Dict{Symbol,Any},
    d::Float64
)
    μ = stats[:mu]
    σ = stats[:sigma]
    exp(-((d - μ)^2) / (2 * σ^2))
end

function predict_membership(
    model::MultiClassModel,
    X::Vector{Float64},
    normalize::Bool=false)

    memberships = Dict{String,Float64}()

    for c in model.classes
        class_freqs = model.class_string_probs[c]
        stats = model.variant_stats[c]
        # d = sum(abs.(X - class_freqs))
        d = sqrt(sum((X - class_freqs) .^ 2))
        memberships[c] = gaussian_membership(stats, d)
    end

    if normalize
        total = sum(values(memberships))
        total > 0 || return memberships
        for variant in keys(memberships)
            memberships[variant] /= total
        end
    end
    return argmax(memberships), memberships
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

function sequence_kmer_distribution(
    regions::Vector{Tuple{Int,Int}},
    seq::Base.CodeUnits,
    kmerset::Vector{String}
)::Vector{UInt64}

    @floop for kmer in kmerset
        local_seq_histogram = zeros(UInt64, length(regions))
        seq_len = length(seq)

        for i in eachindex(regions)
            init_pos, end_pos = regions[i]

            if (end_pos > seq_len)
                end_pos = seq_len
            end

            wndw_buffer = @view seq[init_pos:end_pos]

            if RegionExtraction.occursinKmerBit(codeunits(kmer), wndw_buffer)
                local_seq_histogram[i] += 1
            end
        end

        @reduce(
            kmer_distribution = zeros(UInt64, length(regions)) .+ local_seq_histogram
        )
    end
    return kmer_distribution
end



function predict_raw(
    parameters::Tuple{MultiClassModel,Union{Nothing,String}},
    X::Vector{Float64})::Tuple{String,Dict{String,Float64}}

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

        elseif metric == "rrm"
            # Think in another approach to use RRM, this cross-spectrum approach don't work
            dists[c] = sum(dot(X, class_freqs))

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
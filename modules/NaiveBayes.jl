module NaiveBayes

include("Model.jl")

using FLoops, .Model
export NaiveBayes


struct VariantDataloader
    meta_data::Dict{String,Int}
    byte_seqs::Dict{String,Vector{Base.CodeUnits}}
    wnw_size::Int
end

function getKmerAppearences(
    kmerset::Set{String},
    meta_data::Dict{String,Int},
    byte_seqs::Dict{String,Vector{Base.CodeUnits}},
    wnw_size::Int
)::Vector{Tuple{String,Dict{String,Vector{UInt64}}}}

    dataloader = VariantDataloader(meta_data, byte_seqs, wnw_size)
    get_appearences = Base.Fix2(def_kmer_classes_probs, dataloader)

    distribution = Vector{Tuple{String,Dict{String,Vector{UInt64}}}}(undef, length(kmerset))

    @floop for (i, kmer) in enumerate(collect(kmerset))
        distribution[i] = (kmer, get_appearences(codeunits(kmer)))
    end

    distribution
end

function def_kmer_classes_probs(
    kmer::Base.CodeUnits,
    data::VariantDataloader)::Dict{String,Vector{UInt64}}
    fn_occursin = Base.Fix1(Model.occursinKmerBit, kmer)

    appearance = Dict{String,Vector{UInt64}}()

    for class in keys(data.meta_data)

        seq_windows = maximum(length, data.byte_seqs[class]) - data.wnw_size + 1
        seq_hist = zeros(UInt64, seq_windows)

        for seq in data.byte_seqs[class]
            seq_windows = length(seq) - data.wnw_size + 1

            for initPos in 1:seq_windows
                endPos = initPos + data.wnw_size - 1
                wndw_buffer = @view seq[initPos:endPos]

                if fn_occursin(wndw_buffer)
                    seq_hist[initPos] += 1
                end
            end

        end

        appearance[class] = seq_hist
    end

    appearance
end




end
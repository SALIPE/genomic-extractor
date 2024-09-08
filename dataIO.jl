module DataIO
using FASTX

export DataIO

EIIP_NUCLEOTIDE = Dict{Char,Float64}([
    ('A', 0.1260),
    ('G', 0.0806),
    ('T', 0.1335),
    ('C', 0.1340)])

EIIP_AMINOACID = Dict{Char,Float64}([
    ('L', 0.0000),
    ('I', 0.0000),
    ('N', 0.0036),
    ('G', 0.0050),
    ('V', 0.0057),
    ('E', 0.0058),
    ('P', 0.0198),
    ('H', 0.0242),
    ('K', 0.0371),
    ('A', 0.0373),
    ('Y', 0.0516),
    ('W', 0.0548),
    ('Q', 0.0761),
    ('M', 0.0823),
    ('S', 0.0829),
    ('C', 0.0829),
    ('T', 0.0941),
    ('F', 0.0946),
    ('R', 0.0959),
    ('D', 0.1263)])

function getSequencesFromFastaFile(
    filePath::String
)::Array{FASTX.FASTA.Record}
    sequences::Array{FASTX.FASTA.Record} = []
    for record in open(FASTAReader, filePath)
        push!(sequences, record)
        # push!(sequences, codeunits(sequence(record)))
        # println(identifier(record))
        # println(sequence(record))
        # println(description(record))
    end
    return sequences
end


function getLongestLength(
    filePath::String
)::Int
    sequenceLen::Int = zero(Int)
    for record in open(FASTAReader, filePath)
        if (iszero(sequenceLen))
            sequenceLen = seqsize(record)
        elseif (seqsize(record) > sequenceLen)
            sequenceLen = seqsize(record)
        end

    end
    return seqIndex
end

function getShortestLength(
    filePath::String
)::Int
    sequenceLen::Int = zero(Int)
    for record in open(FASTAReader, filePath)
        @show seqsize(record)
        if (iszero(sequenceLen))
            sequenceLen = seqsize(record)
        elseif (seqsize(record) < sequenceLen)
            sequenceLen = seqsize(record)
        end

    end
    return sequenceLen
end

function sequence2NumericalSerie(
    seqpar::AbstractString
)::Vector{Float64}

    dict = keys(EIIP_NUCLEOTIDE)
    arrSeq = Float64[]
    for c in eachindex(seqpar)
        key = seqpar[c]
        keyin = key ∈ dict
        push!(arrSeq, keyin ? EIIP_NUCLEOTIDE[key] : zero(Float64))
    end
    return arrSeq
end


end
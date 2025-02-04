module DataIO
using FASTX, PyCall

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
)::Tuple{Int,Int}
    seqAmount::Int = zero(Int)
    seqLen::Int = zero(Int)
    for record in open(FASTAReader, filePath)
        seqAmount += 1
        @show seqsize(record)
        if (iszero(seqLen))
            seqLen = seqsize(record)
        elseif (seqsize(record) < seqLen)
            seqLen = seqsize(record)
        end

    end
    return (seqLen, seqAmount)
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

function sequence2NumericalSerie(
    seqpar::AbstractString,
    initIndex::Integer,
    endIndex::Integer
)::Vector{Float64}

    dict = keys(EIIP_NUCLEOTIDE)
    arrSeq = Float64[]
    @inbounds for c in initIndex:endIndex
        key = seqpar[c]
        keyin = key ∈ dict
        push!(arrSeq, keyin ? EIIP_NUCLEOTIDE[key] : zero(Float64))
    end
    return arrSeq
end

function progressBar!(current::Int, total::Int; width::Int=50)
    progress = current / total
    complete_width = Int(round(progress * width))
    incomplete_width = width - complete_width
    bar = "[" * "="^complete_width * "-"^incomplete_width * "]"

    percentage = Int(round(progress * 100))

    print("\r", bar, " ", percentage, "%")

    flush(stdout)
end

function writeFASTASingleChr!(
    originalFilePath::String,
    filename::String,
    ranges::Vector{Tuple{Int,Int}}
)
    for record in open(FASTAReader, originalFilePath)
        name = identifier(record)
        total = length(ranges)
        println("\nExporting file: $name $filename")
        FASTAWriter(open(name * filename, "w")) do writer
            for (current, (initRng, endRng)) in enumerate(ranges)
                DataIO.progressBar!(current, total)
                write(writer, FASTARecord(name * ":" * string(initRng) * "_" * string(endRng),
                    sequence(record)[initRng:endRng]))
            end
        end
    end
end

function writeFASTAS!(
    originalDirPath::String,
    filename::String,
    ranges::Dict{String,Vector{Tuple{Int,Int}}}
)

    chrs = keys(ranges)
    @show chrs
    for file in readdir(originalDirPath)
        println("\nExporting file: $file")
        FASTAWriter(open(file * filename, "w")) do writer
            for record in open(FASTAReader, originalDirPath * "/" * file)
                chr = identifier(record)
                @show chr
                if (chr ∈ chrs)
                    for (initRng, endRng) in ranges[chr]
                        write(writer, FASTARecord(chr * ":" * string(initRng) * "_" * string(endRng),
                            sequence(record)[initRng:endRng]))
                    end
                end
            end
        end
    end

end

function readRegionsFromBed(
    bedFilePath::String
)::Vector{Tuple{String,Int,Int}}

    regioes = Vector{Tuple{String,Int,Int}}()
    open(bedFilePath) do bed
        for linha in readlines(bed)
            if strip(linha) != ""  # Ignora linhas vazias
                campos = split(linha)
                cromossomo::String = campos[1]
                inicio::String = campos[2] # Começo da região (0-indexed)
                fim::String = campos[3]     # Fim da região (1-indexed)
                push!(regioes, (cromossomo, parse(Int, inicio), parse(Int, fim)))
            end
        end

    end
    return regioes
end

function readVectorFromFile(file::String, T::Type)::Vector{T}

    arr = Vector{T}()
    open(file, "r") do fi
        for linha in eachline(fi)
            push!(arr, parse(T, linha))
        end
    end
    return arr
end


function read_pickle_data(file_name)::Vector{String}
    # file_content = read("$variantDirPath/$variant/$(variant)_ExclusiveKmers.txt", String)
    # content_inside_brackets = strip(file_content, ['[', ']'])
    # exclusiveKmers::Vector{String} = strip.(strip.(split(content_inside_brackets, ",")), '\'')
    # data::Vector{String} = pickle.load(open(file_name, "r"), encoding="latin1")

    py"""
import pickle
 
def load_pickle(fpath):
    with open(fpath, "rb") as f:
        data = pickle.load(f)
    return data
"""

    load_pickle = py"load_pickle"(file_name)
    return load_pickle
end


end

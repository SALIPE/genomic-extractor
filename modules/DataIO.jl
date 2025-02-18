module DataIO
using FASTX, Pickle, Serialization, BioSequences

export DataIO

EIIP_NUCLEOTIDE = Dict{Char,Float64}([
    ('A', 0.1260),
    ('G', 0.0806),
    ('T', 0.1335),
    ('C', 0.1340)])

EIIP_AMINOACID = Dict{Char,Float64}([
    ('*', 0.0000),
    ('A', 0.0373),
    ('B', 0.0000),#non-exist
    ('C', 0.0829),
    ('D', 0.1263),
    ('E', 0.0058),
    ('F', 0.0946),
    ('G', 0.0050),
    ('H', 0.0242),
    ('I', 0.0000),
    ('J', 0.0000),#non-exist
    ('K', 0.0371),
    ('L', 0.0000),
    ('M', 0.0823),
    ('N', 0.0036),
    ('O', 0.0000),#non-exist 
    ('P', 0.0198),
    ('Q', 0.0761),
    ('R', 0.0959),
    ('S', 0.0829),
    ('T', 0.0941),
    ('U', 0.0000),#non-exist 
    ('V', 0.0057),
    ('W', 0.0548),
    ('X', 0.0000),#non-exist 
    ('Y', 0.0516),
    ('Z', 0.0000),#non-exist 
])

function padRNA(rna::LongSequence{RNAAlphabet{4}})
    pad_length = (3 - (length(rna) % 3)) % 3
    return rna * LongSequence{RNAAlphabet{4}}(repeat('N', pad_length))
end

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


function sequence2AminNumSerie(
    sequence::String
)::Vector{Float64}

    dna = LongSequence{DNAAlphabet{4}}(sequence)
    rna = convert(LongSequence{RNAAlphabet{4}}, dna)
    amn = BioSequences.translate(padRNA(rna))
    eiip = Vector{Float64}(undef, length(amn))

    @inbounds for (i, c) in enumerate(amn)
        eiip[i] = EIIP_AMINOACID[Char(c)]
    end
    return eiip
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



function writeFASTASingleChr!(
    originalFilePath::String,
    filename::String,
    ranges::Vector{Tuple{Int,Int}}
)
    for record in open(FASTAReader, originalFilePath)
        name = identifier(record)
        println("\nExporting file: $name $filename")
        FASTAWriter(open(name * filename, "w")) do writer
            for (current, (initRng, endRng)) in enumerate(ranges)
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


function read_pickle_data(file_name::AbstractString)
    # file_content = read("$variantDirPath/$variant/$(variant)_ExclusiveKmers.txt", String)
    # content_inside_brackets = strip(file_content, ['[', ']'])
    # exclusiveKmers::Vector{String} = strip.(strip.(split(content_inside_brackets, ",")), '\'')
    # data::Vector{String} = pickle.load(open(file_name, "r"), encoding="latin1")

    load_pickle = Pickle.load(file_name)
    return load_pickle
end

function save_cache(cache_path::String, data)
    try
        open(cache_path, "w") do io
            serialize(io, data)
        end
        @info "Cache saved successfully: $cache_path"
    catch e
        @error "Failed to save cache" exception = (e, catch_backtrace())
    end
end

function load_cache(cache_path::String)
    try
        isfile(cache_path) || return nothing
        open(io -> deserialize(io), cache_path, "r")
    catch e
        @error "Cache loading failed" exception = (e, catch_backtrace())
        nothing
    end
end



end

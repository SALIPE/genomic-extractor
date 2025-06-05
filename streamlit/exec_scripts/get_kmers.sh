#!/bin/bash

if [ $# -ne 3 ]; then
    echo "Uso: $0 <DATA_SOURCE> <WORD_SIZE>"
    echo "Exemplo: $0 ~/datasets/HBV ref 6"
    echo "O script detectará automaticamente os arquivos .fasta no diretório DATA_SOURCE"
    exit 1
fi

source ~/.py-venv/bin/activate
export RUST_BACKTRACE=full

DATA_SOURCE=$1
REF=$2
WORD_SIZE=$3

if [ ! -d "$DATA_SOURCE" ]; then
    echo "Erro: Diretório $DATA_SOURCE não encontrado!"
    exit 1
fi

mkdir -p $DATA_SOURCE/kmers

FASTA_FILES=($(find "$DATA_SOURCE" -maxdepth 1 -name "*.fasta" -type f))

# Verifica se foram encontrados arquivos FASTA
if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "Erro: Nenhum arquivo .fasta encontrado no diretório $DATA_SOURCE"
    exit 1
fi

echo "Arquivos FASTA encontrados: ${#FASTA_FILES[@]}"

for fasta_file in "${FASTA_FILES[@]}"
do
    variant=$(basename "$fasta_file" .fasta)
    
    echo "Processando variante: $variant"
    
    gramep get-only-kmers --rpath $REF \
        --spath "$fasta_file" \
        --save-path $DATA_SOURCE/kmers/ \
        --word $WORD_SIZE \
        --step 1 \
        -d ALL
    
    cp "$fasta_file" $DATA_SOURCE/kmers/$variant
done

echo "Processamento concluído para ${#FASTA_FILES[@]} variantes."
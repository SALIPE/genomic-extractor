#!/bin/bash
source ~/.py-venv/bin/activate

export RUST_BACKTRACE=full

DENGUE=~/Desktop/datasets/dengue
HBV=~/Desktop/datasets/HBV/data

GREAC=~/Desktop/genomic-extractor/scripts/local/benchmark.sh
BALANCEDATASET=~/Desktop/Fasta-splitter/FastaSplitter

REF_HIV=../genomic-extractor/comparison_scripts/castor_hiv_data/hiv1_refseq.fasta
REF_HBV=~/Desktop/datasets/HBV/refseq.fasta
REF_DENV=~/Desktop/datasets/denv/refseq.fasta

if [ $# -lt 2 ]; then
    echo "❌ Erro: Argumentos insuficientes"
    echo "Uso: $0 <GROUPNAME> <WINDOW>"
    exit 1
fi

GROUPNAME=$1
WINDOW=$2

echo "📋 Parâmetros recebidos:"
echo "   - GROUPNAME: $GROUPNAME"
echo "   - WINDOW: $WINDOW"

case $GROUPNAME in
    denv)
        SOURCE=$DENGUE
        echo "✅ Dataset DENGUE selecionado: $SOURCE"
        ;;
    hbv)
        SOURCE=$HBV
        echo "✅ Dataset HBV selecionado: $SOURCE"
        ;;
    *)
        echo "❌ Erro: GROUPNAME deve ser 'denv' ou 'hbv'"
        exit 1
        ;;
esac


function get_kmers_denv() {
    
    for variant in type1 type2 type3 type4; do
        
        gramep get-only-kmers \
            --rpath $REF_DENV \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word 6 \
            --step 1
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

function get_kmers_hbv() {
    
    for variant in A B C D E F; do
        gramep get-only-kmers \
            --rpath $REF_HBV \
            --spath $SOURCE/train/$variant.fasta \
            --save-path $SOURCE/train/kmers/ \
            --word 6 \
            --step 1
        
        mv $SOURCE/train/$variant.fasta kmers/$variant/$variant.fasta
    done
}

TRAIN=$SOURCE/train/kmers
TESTDIR=$SOURCE/test
METRIC=manhattan

echo "📁 Caminhos configurados:"
echo "   - TRAIN: $TRAIN"
echo "   - TESTDIR: $TESTDIR"
echo "   - METRIC: $METRIC"

echo "🔍 Verificando existência dos diretórios..."
if [ ! -d "$SOURCE" ]; then
    echo "❌ Erro: Diretório SOURCE não existe: $SOURCE"
    exit 1
fi

if [ ! -d "$SOURCE/train" ]; then
    echo "❌ Erro: Diretório de treino não existe: $SOURCE/train"
    exit 1
fi

if [ ! -d "$TESTDIR" ]; then
    echo "❌ Erro: Diretório de teste não existe: $TESTDIR"
    exit 1
fi

echo "✅ Todos os diretórios verificados"

for i in {1..100}; do
    echo "Iteração $i de 100"
    
    $BALANCEDATASET/test.sh $SOURCE

    cd $SOURCE/train
    mkdir -p kmers

     case $GROUPNAME in
        denv)
            get_kmers_denv
            ;;
        hbv)
            get_kmers_hbv
            ;;
        *)
            echo "❌ Erro: GROUPNAME inválido: $GROUPNAME"
            exit 1
            ;;
    esac
    
    $GREAC $TRAIN $TESTDIR $GROUPNAME $WINDOW $METRIC 
    
done


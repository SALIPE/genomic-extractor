import json
import os
import shutil
import subprocess
import zipfile

import pandas as pd
import plotly.express as px
import requests
from Bio import SeqIO

import streamlit as st

# Configuração da página
st.set_page_config(
    page_title="GREAC-UI",
    page_icon="🧬",
    layout="wide"
)

def list_directory_tree(directory: str, prefix: str = "") -> str:
    tree_str = ""
    for root, dirs, files in os.walk(directory):
        level = root.replace(directory, "").count(os.sep)
        indent = " " * 4 * level
        tree_str += f"{prefix}{indent}📁 {os.path.basename(root)}/\n"
        sub_indent = " " * 4 * (level + 1)
        for f in files:
            tree_str += f"{prefix}{sub_indent}📄 {f}\n"
    return tree_str

INPUT_DIR = "/app/data/input"
OUTPUT_DIR = "./data/output"
JULIA_CACHE_DIR = "/julia-cache"
HOMEDATA = "/home/salipe/"

def extract_uploaded_files(uploaded_files):
    """Extrai arquivos uploaded para o diretório de input"""
    
    for uploaded_file in uploaded_files:
        file_path = os.path.join(INPUT_DIR, uploaded_file.name)
        
        if uploaded_file.name.endswith('.zip'):
            # Extrai arquivo ZIP
            with open(file_path, "wb") as f:
                f.write(uploaded_file.getbuffer())
            
            with zipfile.ZipFile(file_path, 'r') as zip_ref:
                zip_ref.extractall(INPUT_DIR)
            
            os.remove(file_path)  # Remove o arquivo ZIP após extração
        else:
            raise Exception("Model file not structure")

def analyze_input_structure():
    """Analisa a estrutura dos arquivos de input"""
    variants = {}
    
    if not os.path.exists(INPUT_DIR):
        return variants
    
    for root, dirs, files in os.walk(INPUT_DIR):
        for file in files:
            if file.endswith(('.fasta', '.fa', '.fas')):
                variant_name = os.path.basename(root) if root != INPUT_DIR else "root"
                file_path = os.path.join(root, file)
                
                if variant_name not in variants:
                    variants[variant_name] = []
                
                # Conta sequências no arquivo FASTA
                try:
                    seq_count = len(list(SeqIO.parse(file_path, "fasta")))
                    variants[variant_name].append({
                        'file': file,
                        'path': file_path,
                        'sequences': seq_count
                    })
                except:
                    variants[variant_name].append({
                        'file': file,
                        'path': file_path,
                        'sequences': 0
                    })
    
    return variants



def main():
    st.title("🧬 GREAC-UI")
    st.markdown("---")

    
    st.header("⚙️ Configurações")

    col1, col2 = st.columns(2)
    with col1:
        group_name = st.text_input("Nome do grupo")
        window_size = st.number_input("Comprimento mínimo da sequência", step=0.001, format="%.3f")
        metric = st.selectbox("Métricas", ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        cache = st.checkbox("Usar cache", value=True)
    with col2:
        pg_options = ["GREAC", "gramep", "FastaSplitter"]
        program = st.selectbox("Programa",pg_options)

        if program == pg_options[1]:
            reference_path = st.text_input("📂 Caminho da Referência", value=f"{HOMEDATA}/Desktop/datasets/denv/refseq.fasta")
            k = st.number_input("K-mer size", step=1)

    st.header("📁 Seleção de Diretórios Locais")

    col1, col2, col3 = st.columns(3)
    
    with col1:
        train_dir = st.text_input("📂 Caminho dos Dados de Treino", value=f"{HOMEDATA}/Desktop/datasets/dengue/train/kmers")
    with col2:
        test_dir = st.text_input("📂 Caminho dos Dados de Teste", value=f"{HOMEDATA}/Desktop/datasets/dengue/test")
    with col3:
        feature_dir = st.text_input("📂 Caminho dos Dados para Extração de Features", value=f"{HOMEDATA}/Desktop/datasets/dengue")

    # Verificação dos diretórios
    if program == pg_options[0]:
        dir_check = all([os.path.isdir(train_dir), os.path.isdir(test_dir)])
    else: 
        dir_check = all([os.path.isdir(feature_dir)])

    if not dir_check:
        st.warning("⚠️ Um ou mais diretórios não existem. Verifique os caminhos.")
    else:
        st.success("✅ Diretórios encontrados.")
        st.header("🚀 Processamento")

        # Campos obrigatórios validados
        if (not group_name or not window_size) and program == pg_options[0]:
            st.error("⚠️ Preencha todos as configurações antes de iniciar o processamento.")
            return

        if program == pg_options[0]:
            process_script = f'../scripts/local/benchmark.sh'
        elif program == pg_options[1]:
            process_script = './scripts/get_kmers.sh'
        elif program == pg_options[2]:
            process_script = f'{HOMEDATA}/Desktop/Fasta-splitter/FastaSplitter/balance.sh'

        if not os.path.exists(process_script):
            st.error("❌ Script de processamento não encontrado.")
            return

        col1, col2 = st.columns(2)
            
        with col1:
            if st.button("▶️ Iniciar Processamento", type="primary"):
                st.info("📡 Iniciando o script de processamento...")
                if program == pg_options[0]:
                    st.text(f"Grupo: {group_name} | Window: {window_size} | Métrica: {metric}")

                greac_cmd = [
                    process_script,
                    train_dir,
                    test_dir,
                    group_name,
                    str(window_size),
                    metric,
                    "--no-cache" if not cache else ""
                ]
                balance = [
                    process_script,
                    feature_dir
                ]

                get_kmers = [
                    process_script,
                    train_dir,
                    reference_path,
                    k
                ]

                if program == pg_options[0]:
                    process_cmd = greac_cmd
                elif program == pg_options[1]:
                    process_cmd = get_kmers
                else: 
                    process_cmd = balance

                with st.spinner("🔄 Processando..."):
                    process = subprocess.Popen(
                        process_cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        universal_newlines=True
                    )

                    output_placeholder = st.empty()
                    output_text = ""

                    for line in process.stdout:
                        output_text += line
                        output_placeholder.code(output_text, language='bash')

                    process.wait()
                    st.success("✅ Processamento concluído.")
            
        with col2:
            if st.button("📥 Baixar Resultados"):
                        if os.path.exists(OUTPUT_DIR) and os.listdir(OUTPUT_DIR):
                            # Cria ZIP com resultados
                            zip_path = "/tmp/results.zip"
                            with zipfile.ZipFile(zip_path, 'w') as zipf:
                                for root, dirs, files in os.walk(OUTPUT_DIR):
                                    for file in files:
                                        file_path = os.path.join(root, file)
                                        arcname = os.path.relpath(file_path, OUTPUT_DIR)
                                        zipf.write(file_path, arcname)
                            
                            with open(zip_path, "rb") as f:
                                st.download_button(
                                    label="💾 Download ZIP",
                                    data=f.read(),
                                    file_name="fasta_results.zip",
                                    mime="application/zip"
                                )
                        else:
                            st.warning("⚠️ Nenhum resultado encontrado")



if __name__ == "__main__":
    main()
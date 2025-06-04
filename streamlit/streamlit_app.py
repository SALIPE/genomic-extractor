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
OUTPUT_DIR = "/app/data/output"
JULIA_CACHE_DIR = "/julia-cache"
N8N_WEBHOOK_URL = "http://n8n:5678/webhook/process-fasta"

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

def trigger_n8n_workflow(processing_params):
    """Dispara o workflow do N8N"""
    try:
        response = requests.post(
            N8N_WEBHOOK_URL,
            json=processing_params,
            timeout=None
        )
        return response.status_code == 200, response.json() if response.content else {}
    except Exception as e:
        return False, str(e)


def main():
    st.title("🧬 GREAC-UI")
    st.markdown("---")

    with st.sidebar:
        st.header("⚙️ Configurações")

        group_name = st.text_input("Nome do grupo")
        window_size = st.number_input("Comprimento mínimo da sequência", step=0.005, format="%.3f")
        metric = st.selectbox("Métricas", ["manhattan", "euclidian", "chisquared", "mahalanobis", "kld"])
        cache = st.checkbox("Usar cache", value=True)

    st.header("📁 Seleção de Diretórios Locais")

    train_dir = st.text_input("📂 Caminho dos Dados de Treino", value="./data/train")
    test_dir = st.text_input("📂 Caminho dos Dados de Teste", value="./data/test")
    feature_dir = st.text_input("📂 Caminho dos Dados para Extração de Features", value="./data/feature")

    # Verificação dos diretórios
    dir_check = all([os.path.isdir(train_dir), os.path.isdir(test_dir), os.path.isdir(feature_dir)])

    if not dir_check:
        st.warning("⚠️ Um ou mais diretórios não existem. Verifique os caminhos.")
    else:
        st.success("✅ Diretórios encontrados.")

        st.header("🚀 Processamento")

        # Campos obrigatórios validados
        if not group_name or not window_size:
            st.error("⚠️ Preencha todos os campos do lado esquerdo antes de iniciar o processamento.")
            return

        process_script = '../scripts/local/benchmark.sh'

        if not os.path.exists(process_script):
            st.error("❌ Script de processamento não encontrado.")
            return

        if st.button("▶️ Iniciar Processamento", type="primary"):
            st.info("📡 Iniciando o script de processamento...")
            st.text(f"Grupo: {group_name} | Window: {window_size} | Métrica: {metric}")

            
            cmd = [
                process_script,
                train_dir,
                test_dir,
                group_name,
                str(window_size),
                metric,
                "--no-cache" if not cache else ""
            ]

            with st.spinner("🔄 Processando..."):
                process = subprocess.Popen(
                    cmd,
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



if __name__ == "__main__":
    main()
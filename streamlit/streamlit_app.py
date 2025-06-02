import json
import os
import shutil
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

# Diretórios
INPUT_DIR = "/app/data/input"
OUTPUT_DIR = "/app/data/output"
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

def check_processing_status():
    """Verifica status do processamento"""
    if os.path.exists(os.path.join(OUTPUT_DIR, "processing_complete.json")):
        with open(os.path.join(OUTPUT_DIR, "processing_complete.json"), 'r') as f:
            return json.load(f)
    return None

def main():
    st.title("🧬 FASTA Variant Processor")
    st.markdown("---")
    
    # Sidebar para configurações
    with st.sidebar:
        st.header("⚙️ Configurações")
        
        processing_mode = st.selectbox(
            "Modo de Processamento",
            ["benchmark", "fit-parameters"]
        )

        group_name = st.text_input(
            "Nome do grupo"
        )

        window_size = st.number_input(
            "Comprimento mínimo da sequência",
            step=.05,format="%.3f"
        )
           
    # Upload de arquivos
    st.header("📁 Upload de Arquivos")
    uploaded_files = st.file_uploader(
        "Selecione arquivos ZIP contendo estrutura de variantes",
        accept_multiple_files=True,
        type=['zip']
    )
    
    if uploaded_files:
        with st.spinner("Extraindo arquivos..."):
            extract_uploaded_files(uploaded_files)
        st.success(f"✅ {len(uploaded_files)} arquivo(s) processado(s)")
        
        # Análise da estrutura
        st.header("🔍 Estrutura dos Dados")
        variants = analyze_input_structure()
        
        if variants:
            col1, col2 = st.columns(2)
            
            with col1:
                st.metric("Variantes Encontradas", len(variants))
                
                # Tabela de variantes
                variant_data = []
                for variant, files in variants.items():
                    total_sequences = sum(f['sequences'] for f in files)
                    variant_data.append({
                        'Variante': variant,
                        'Arquivos': len(files),
                        'Sequências': total_sequences
                    })
                
                df_variants = pd.DataFrame(variant_data)
                st.dataframe(df_variants, use_container_width=True)
            
            with col2:
                # Gráfico de distribuição
                if len(variant_data) > 0:
                    fig = px.bar(
                        df_variants,
                        x='Variante',
                        y='Sequências',
                        title='Distribuição de Sequências por Variante'
                    )
                    st.plotly_chart(fig, use_container_width=True)
            
            # Detalhes por variante
            st.subheader("📋 Detalhes por Variante")
            for variant, files in variants.items():
                with st.expander(f"Variante: {variant} ({len(files)} arquivos)"):
                    for file_info in files:
                        st.write(f"- **{file_info['file']}**: {file_info['sequences']} sequências")
            
            # Botão de processamento
            st.header("🚀 Processamento")
            
            processing_params = {
                "group_name": group_name,
                "window_size":window_size,
                "train_dir": f'{INPUT_DIR}/train/kmers',
                "test_dir": f'{INPUT_DIR}/test'
            }
            
            col1, col2, col3 = st.columns(3)
            
            with col1:
                if st.button("▶️ Iniciar Processamento", type="primary"):
                    with st.spinner("Iniciando processamento..."):
                        success, result = trigger_n8n_workflow(processing_params)
                        
                        if success:
                            st.success("✅ Processamento iniciado com sucesso!")
                            st.json(result)
                        else:
                            st.error(f"❌ Erro ao iniciar processamento: {result}")
            
            with col2:
                if st.button("🔄 Verificar Status"):
                    status = check_processing_status()
                    if status:
                        st.success("✅ Processamento concluído!")
                        st.json(status)
                    else:
                        st.info("⏳ Processamento em andamento ou não iniciado")
            
            with col3:
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
    
    # Status do sistema
    st.header("📊 Status do Sistema")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        # Verifica conexão com N8N
        try:
            response = requests.get("http://n8n:5678/healthz", timeout=5)
            n8n_status = "🟢 Online" if response.status_code == 200 else "🔴 Offline"
        except:
            n8n_status = "🔴 Offline"
        st.metric("N8N Status", n8n_status)
    
    with col2:
        input_files = len([f for f in os.listdir(INPUT_DIR) if f.endswith(('.fasta', '.fa', '.fas'))]) if os.path.exists(INPUT_DIR) else 0
        st.metric("Arquivos Input", input_files)
    
    with col3:
        output_files = len(os.listdir(OUTPUT_DIR)) if os.path.exists(OUTPUT_DIR) else 0
        st.metric("Arquivos Output", output_files)

if __name__ == "__main__":
    main()
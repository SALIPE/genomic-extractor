FROM python:3.9-slim

WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first for better caching
COPY streamlit_requirements.txt .
RUN pip install --no-cache-dir -r streamlit_requirements.txt

# Copy application code
COPY . .

# Create necessary directories
RUN mkdir -p /app/data/input /app/data/output

EXPOSE 8501

CMD ["streamlit", "run", "streamlit_app.py", "--server.address", "0.0.0.0"]
FROM python:3.12-slim

# Dependências essenciais compatíveis com Debian Trixie
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libffi-dev \
    libssl-dev \
    curl \
    libgl1 \
    libglib2.0-0 \
    libc6 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app.py .
COPY degs/ ./degs/
COPY go/ ./go/
COPY images/ ./images/
COPY lib/ ./lib/
COPY kegg/ ./kegg/
COPY supporting/ ./supporting/
COPY analysis/ ./analysis/

EXPOSE 8505

ENV STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_SERVER_ENABLECORS=false \
    STREAMLIT_SERVER_ENABLE_XSRF_PROTECTION=false \
    STREAMLIT_SERVER_ADDRESS=0.0.0.0 \
    STREAMLIT_SERVER_PORT=8501 \
    STREAMLIT_BROWSER_GATHERUSAGESTATS=false

CMD ["streamlit", "run", "app.py"]
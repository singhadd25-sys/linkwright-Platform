FROM python:3.11-slim
ENV PYTHONUNBUFFERED=1 LINKWRIGHT_DATA=/data
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends libexpat1 libxrender1 libxext6 && rm -rf /var/lib/apt/lists/*
COPY backend/requirements.txt backend/requirements.txt
RUN pip install --no-cache-dir -r backend/requirements.txt
COPY backend backend
COPY frontend frontend
RUN useradd -m app && mkdir /data && chown app /data
USER app
WORKDIR /app/backend
EXPOSE 8000
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]

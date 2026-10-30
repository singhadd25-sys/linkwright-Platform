# Linkwright — Computational Molecular Research

This repository contains the public GitHub Pages landing site plus the Linkwright v0.4.2 FastAPI/RDKit research application.

- Public site: `https://singhadd25-sys.github.io/linkwright-Platform/`
- Research UI: `https://singhadd25-sys.github.io/linkwright-Platform/app/`
- Backend: FastAPI + RDKit + AutoDock Vina integration + GROMACS package preparation
- Authentication: Supabase Auth; API routes require a valid bearer token

## Deployment
See [DEPLOYMENT.md](DEPLOYMENT.md). The static UI can be served by GitHub Pages. Computational endpoints require the Docker backend to be deployed separately.

## Local development
1. Copy `.env.example` to `.env` and set the Supabase URL and publishable key.
2. Run `docker compose up -d --build`.
3. Open `http://localhost:8000` for the local FastAPI-served frontend.

Never commit `.env`, service-role keys, private research data, or user job outputs.

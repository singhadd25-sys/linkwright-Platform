# Linkwright deployment and authentication

## Architecture
- GitHub Pages serves the public marketing site and authenticated research UI at `/app/`.
- Supabase Auth manages researcher sign-in.
- The FastAPI/RDKit Docker API is a separate HTTPS service (Render blueprint included).
- The API validates each Supabase bearer token against Supabase Auth. Job status, result and download routes enforce owner IDs.
- The browser contains only the Supabase publishable key. Never expose a Supabase secret/service-role key.

## Deploy the backend (one-time external step)
1. Push this repository to GitHub.
2. Open Render and create a new Blueprint from `singhadd25-sys/linkwright-Platform`; select `render.yaml`.
3. In the Render service environment, set `SUPABASE_URL` to `https://ppxrlrevcitafcksajpi.supabase.co` and `SUPABASE_PUBLISHABLE_KEY` to the project's publishable key.
4. Deploy. Check `https://YOUR-RENDER-SERVICE.onrender.com/api/health` returns `{"ok":true,...}`.
5. Replace the placeholder `LINKWRIGHT_API_URL` in `app/config.js` with that HTTPS service URL, commit and push.

The Render free web service may sleep when idle and its filesystem is ephemeral. This prototype stores job files on the service filesystem. Do not use it for valuable or confidential research data until private persistent object storage and a retention/backup plan are configured. The API is intentionally not marked production-ready for multi-user long-running job execution.

## Supabase Auth URLs
Set Site URL to `https://singhadd25-sys.github.io/linkwright-Platform/` and add these exact Redirect URLs:
- `https://singhadd25-sys.github.io/linkwright-Platform/`
- `https://singhadd25-sys.github.io/linkwright-Platform/app/`

## Local Docker
Copy `.env.example` to `.env`, set the Supabase URL and publishable key, then run `docker compose up -d --build`. Never commit `.env`.

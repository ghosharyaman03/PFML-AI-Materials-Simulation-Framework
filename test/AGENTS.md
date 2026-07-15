# AGENTS.md

## What this repo is

Docker-based reproducible dev environment for scientific computing + AI agent orchestration via AutoGen + Ollama.

## Key files

| File | Purpose |
|------|---------|
| `Dockerfile` | Ubuntu 22.04 image with Python, Rust, scientific C/C++ libs (OpenBLAS, FFTW, Eigen, Boost, VTK), PyTorch (CPU), TensorFlow, JAX, OpenCode CLI |
| `docker-compose.yml` | Runs container with `/app` (rw) and `../` as `/shared` (ro); exposes `host.docker.internal` for host Ollama |
| `setup.sh` | Build image (`docker build -t msbase .`) then run interactively with 6G memory limit, mounts, and pids-limit=500 |
| `run_agents.py` | Standalone AutoGen script: two agents (user_proxy + orchestrator) talking to `host.docker.internal:11434/v1` (gemma2:2b), outputs to `/app/test` |
| `test/` | Empty dir; used as `work_dir` by `run_agents.py` |

## Commands

```bash
# Build and run container interactively
bash setup.sh

# Or via compose (stays running in background)
docker compose up -d

# Run AutoGen agent loop inside container
python run_agents.py
```

## Architecture notes

- `run_agents.py` connects to host's Ollama via `host.docker.internal:11434` — requires Ollama running on host with `gemma2:2b` pulled
- No build/test/lint tooling (pytest is installed in image but no tests exist yet)
- `/shared` on host is mounted read-only for data access; all work happens in `/app`
- The `test/` directory is the code execution sandbox for agents
- Image includes OpenCode CLI (`opencode`) for future agent usage
- The "Google Antigravity" entries in the Dockerfile are a reference to XKCD #353 and not a real service

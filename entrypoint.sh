#!/bin/bash
# entrypoint.sh
# ============================================================
# 1. Starts the Ollama server in the background.
# 2. Verifies your source was actually mounted to /app (fails fast
#    with a clear message instead of a confusing Python import
#    error later if you forgot -v $(pwd):/app).
# 3. Bakes/refreshes the pfml-parser model from the mounted
#    build_modelfile.py + shared_core.py + Example_Systems, so it
#    always reflects whatever's currently on disk -- including any
#    edits made since the image was built.
# 4. Ensures the DATA/ output directory exists.
# 5. Hands off to whatever command the container was started with.
# ============================================================
set -euo pipefail

ollama serve > /var/log/ollama.log 2>&1 &

# Bounded wait for a local, fast, deterministic startup -- not a
# flaky network operation, so a simple timeout-bounded poll is the
# right tool here, not a hand-rolled retry-and-count loop.
if ! timeout 30 bash -c 'until curl -sf http://localhost:11434/ >/dev/null 2>&1; do sleep 0.5; done'; then
    echo "entrypoint.sh: Ollama server did not become ready within 30s" >&2
    echo "entrypoint.sh: see /var/log/ollama.log for details" >&2
    exit 1
fi

# --- Verify the source mount actually landed ---
if [ ! -f /app/build_modelfile.py ] || [ ! -f /app/shared_core.py ]; then
    echo "entrypoint.sh: /app/build_modelfile.py or /app/shared_core.py not found." >&2
    echo "entrypoint.sh: did you forget to mount your source? e.g.:" >&2
    echo "entrypoint.sh:   docker run -it -v \$(pwd):/app msbase" >&2
    exit 1
fi

# --- Bake pfml-parser from the mounted source ---
echo "entrypoint.sh: baking pfml-parser from mounted Example_Systems..."
python3 /app/build_modelfile.py --base-model "${OLLAMA_BASE_MODEL}" --out-model pfml-parser

# --- Ensure output directory exists (mount may not include it yet) ---
mkdir -p "${MICROSIM_DIR:-/app/microsim_gp}/DATA"

exec "$@"

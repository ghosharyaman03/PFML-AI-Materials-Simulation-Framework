#!/bin/bash
# Starts the Ollama server in the background, waits until it's responding,
# then (re)builds the pfml-parser custom model from the CURRENT contents
# of Example_Systems/ before handing off to whatever command the
# container was started with (default: bash).
#
# This runs at container START, not image build time, on purpose: it's
# what lets pfml-parser stay in sync with presets added or edited via a
# live-mounted volume after the image was built, instead of silently
# working off a stale system list baked in at `docker build` time.
set -e

ollama serve > /var/log/ollama.log 2>&1 &

for i in $(seq 1 30); do
    if curl -sf http://localhost:11434/ >/dev/null 2>&1; then
        break
    fi
    sleep 1
done

# Build/refresh the custom model (base model must already be pulled --
# that happens once at image build time in the Dockerfile, since it's
# large and doesn't depend on live data). If this fails for any reason
# (e.g. base model missing), don't block the container from starting --
# stage1_translator.py falls back to regex-only parsing if Ollama/the
# model isn't reachable.
if [ -f /app/build_modelfile.py ]; then
    python3 /app/build_modelfile.py \
        --base-model "${OLLAMA_BASE_MODEL:-gemma3:4b}" \
        --out-model "${OLLAMA_MODEL:-pfml-parser}" \
        || echo "  [entrypoint] Warning: pfml-parser build failed -- falling back to base model/regex parsing."
fi

exec "$@"

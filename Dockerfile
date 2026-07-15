FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC
ENV PYTHONUNBUFFERED=1

SHELL ["/bin/bash", "-c"]

RUN echo 'Acquire::Retries "5";' > /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::http::Timeout "30";' >> /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::https::Timeout "30";' >> /etc/apt/apt.conf.d/80-retries
# ============================================================
# System + Build Tools
# ============================================================
RUN for i in 1 2 3 4 5; do \
      apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ gfortran \
        make pkg-config cmake \
        git curl wget \
        ca-certificates \
        locales \
      && break || { echo "apt attempt $i failed, retrying..."; sleep 5; }; \
    done && rm -rf /var/lib/apt/lists/*
# ============================================================
# Scientific Libraries (HDF5, GSL — required by microsim_gp)
# ============================================================
RUN for i in 1 2 3 4 5; do \
      apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-openmpi-dev hdf5-tools libgsl-dev libgsl27 \
      && break || { echo "apt attempt $i failed, retrying..."; sleep 5; }; \
    done && ldconfig && rm -rf /var/lib/apt/lists/*
# ============================================================
# MPI
# ============================================================
RUN for i in 1 2 3 4 5; do \
      apt-get update && apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi-dev \
      && break || { echo "apt attempt $i failed, retrying..."; sleep 5; }; \
    done && rm -rf /var/lib/apt/lists/*


# ============================================================
# Ollama (local LLM runtime for the NL translator stage)
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends zstd && rm -rf /var/lib/apt/lists/*

RUN for i in 1 2 3 4 5; do \
      curl -fL https://ollama.com/download/ollama-linux-amd64.tar.zst -o /tmp/ollama.tar.zst \
      && tar --zstd -C /usr -xf /tmp/ollama.tar.zst \
      && rm -f /tmp/ollama.tar.zst \
      && break || { echo "ollama install attempt $i failed, retrying..."; sleep 5; }; \
    done && \
    ollama --version


# Pull the model at build time so it's baked into the image and doesn't
# need a network call at container start. This starts a temporary
# ollama server, waits for it to be ready, pulls the model, then stops
# the server again -- the actual runtime server is started by the
# entrypoint script below.
ENV OLLAMA_MODEL=gemma3:4b
RUN ollama serve & \
    OLLAMA_PID=$! && \
    for i in $(seq 1 30); do \
      curl -sf http://localhost:11434/ >/dev/null 2>&1 && break || sleep 1; \
    done && \
    timeout 3600 ollama pull ${OLLAMA_MODEL} && \
    kill $OLLAMA_PID


# ============================================================
# Python
# ============================================================
RUN for i in 1 2 3 4 5; do \
      apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-dev python-is-python3 \
      && break || { echo "apt attempt $i failed, retrying..."; sleep 5; }; \
    done && rm -rf /var/lib/apt/lists/*

RUN pip3 install --no-cache-dir \
    numpy \
    pycalphad

# ============================================================
# Copy project files
# ============================================================
WORKDIR /app
COPY microsim_gp/ /app/microsim_gp/
COPY shared_core.py /app/shared_core.py
COPY stage1_translator.py /app/stage1_translator.py
COPY stage2_builder.py /app/stage2_builder.py
COPY stage3_validator.py /app/stage3_validator.py
COPY run_agents.py /app/run_agents.py
COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh
COPY build_modelfile.py /app/build_modelfile.py

# Bake the valid-systems/parameter knowledge into a custom Ollama model
# ("pfml-parser") on top of the base model, so stage1_translator.py sends
# only the user's raw prompt per call instead of the full systems list
# every time. To upgrade later (e.g. to a 4B/8B model), just change
# OLLAMA_BASE_MODEL below -- no code changes needed elsewhere.
ENV OLLAMA_BASE_MODEL=gemma3:4b
ENV OLLAMA_MODEL=pfml-parser

# --- Verify the source actually made it into the image ---
RUN echo "===== Contents of /app/microsim_gp =====" \
    && ls -la /app/microsim_gp \
    && find /app/microsim_gp -maxdepth 2 -iname "*.c" \
    && test -f /app/microsim_gp/microsim_gp.c \
        || (echo "FATAL: microsim_gp.c not found under /app/microsim_gp/. Check your host folder layout." && exit 1)

WORKDIR /app/microsim_gp

# --- Compile, and fail the build hard if it doesn't work ---
RUN HDF5_CFLAGS=$(pkg-config --cflags hdf5-openmpi 2>/dev/null || echo "-I/usr/include/hdf5/openmpi") && \
    HDF5_LIBS=$(pkg-config --libs hdf5-openmpi 2>/dev/null || echo "-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5") && \
    echo "===== Compiling: mpicc microsim_gp.c $HDF5_CFLAGS $HDF5_LIBS =====" && \
    mpicc -o /usr/local/bin/microsim_gp microsim_gp.c \
        -Ifunctions -Isolverloop \
        $HDF5_CFLAGS $HDF5_LIBS \
        -lm -lgsl -lgslcblas \
    || (echo "FATAL: mpicc compilation failed. See errors above." && exit 1)

# --- Verify the binary actually exists before declaring victory ---
RUN test -x /usr/local/bin/microsim_gp \
        || (echo "FATAL: binary missing despite mpicc reporting success." && exit 1) \
    && ln -sf /usr/local/bin/microsim_gp /app/microsim_gp/microsim_gp \
    && ls -la /usr/local/bin/microsim_gp /app/microsim_gp/microsim_gp \
    && echo "microsim_gp compiled successfully"

# ============================================================
# Create DATA output directory
# ============================================================
RUN mkdir -p /app/microsim_gp/DATA

# ============================================================
# Environment
# ============================================================
ENV MICROSIM_DIR=/app/microsim_gp
ENV MPI_NP=4
ENV MPI_GRID="2 2"

WORKDIR /app
ENTRYPOINT ["entrypoint.sh"]
CMD ["bash"]

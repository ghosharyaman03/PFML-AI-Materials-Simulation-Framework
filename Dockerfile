# syntax=docker/dockerfile:1
#
# PFML / msbase -- production Dockerfile (mount-based dev variant)
# ============================================================
# This image is a reusable ENVIRONMENT + TOOLCHAIN only. Application
# source (microsim_gp/, shared_core.py, stage*.py, run_agents.py,
# build_modelfile.py) is NOT baked in -- it's bind-mounted at
# /app at `docker run` time instead, so edits on the host take
# effect immediately with no rebuild. Compiling microsim_gp.c is
# likewise a manual, on-demand step against the mounted source, not
# something this image does for you.
#
# Because application code only exists once the volume is mounted,
# anything that depends on it (baking the pfml-parser Ollama model
# from build_modelfile.py + shared_core.py) has to happen at
# container START, not image BUILD -- see entrypoint.sh. Only the
# raw base-model download stays at build time, since that's a large,
# slow pull with no dependency on your source code.
#
# Ordering principle: every layer is installed in the order its
# dependencies require it, so nothing is ever installed before
# something it needs is already present.
#
# Retry philosophy: no hand-rolled "for i in 1 2 3 4 5; do ... done"
# loops anywhere in this file. apt and curl both have native retry
# support -- we configure and use that instead, once, correctly.
#
# Configure the base model and MPI topology at build time with
# --build-arg instead of editing this file:
#   docker build --build-arg OLLAMA_BASE_MODEL=gemma3:4b -t msbase .
#
# Run with your source mounted in, e.g.:
#   docker run -it -v $(pwd):/app --gpus all msbase
# ============================================================

FROM ubuntu:22.04

ARG OLLAMA_BASE_MODEL=gemma2:2b
ARG MPI_NP=4
ARG MPI_GRID="2 2"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Every RUN below fails fast and loud on any error -- no step can
# silently succeed despite an internal command failing.
SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# ============================================================
# 1. APT reliability configuration
#    apt's own retry/timeout handling, configured once, replaces
#    every hand-rolled "for i in 1 2 3 4 5" apt loop.
# ============================================================
RUN echo 'Acquire::Retries "5";'              > /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::http::Timeout "30";'       >> /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::https::Timeout "30";'      >> /etc/apt/apt.conf.d/80-retries

# ============================================================
# 2. Core system utilities
#    Needed by every step after this one: curl to fetch Ollama,
#    ca-certificates so HTTPS fetches actually validate.
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# 3. Compiler / build toolchain
#    You compile microsim_gp.c yourself against the mounted source,
#    but the toolchain itself still belongs in the image -- it's
#    part of the reusable environment, not the application code.
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential gcc g++ gfortran \
        make pkg-config cmake \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# 4. MPI runtime + headers
#    Installed before HDF5, because the HDF5 variant we need
#    (libhdf5-openmpi-dev) is built against it.
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        openmpi-bin libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# 5. Numerical libraries (depend on MPI already being present)
#    HDF5 (MPI-aware variant), GSL (spline interpolation, used by
#    Function_F=4's CALPHAD-backed free energy), BLAS/LAPACK.
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-openmpi-dev hdf5-tools \
        libgsl-dev libgsl27 \
        libopenblas-dev liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# ============================================================
# 6. Python interpreter + package manager
#    Needed before any pip install, and before build_modelfile.py /
#    the pipeline stages can run later (once mounted in).
# ============================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 python3-pip python3-dev python-is-python3 \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --upgrade pip

# ============================================================
# 7. Python packages
#    Installed after the BLAS/LAPACK/HDF5 system libraries above so
#    numpy/scipy/h5py link against them rather than pulling in
#    redundant bundled copies.
# ============================================================
RUN pip3 install --no-cache-dir \
        numpy \
        pandas \
        scipy \
        pycalphad \
        h5py \
        matplotlib

# ============================================================
# 8. Ollama runtime + base model
#    The binary install and the raw model pull both belong at build
#    time: neither depends on your application source, and the
#    model weights are a large, slow download you don't want
#    repeated on every container start. curl's own --retry flags
#    replace the hand-rolled download loop.
# ============================================================
RUN curl -fL --retry 50 --retry-delay 5 --retry-all-errors --retry-connrefused \
        -C - --connect-timeout 30 \
        https://ollama.com/download/ollama-linux-amd64.tar.zst \
        -o /tmp/ollama.tar.zst \
    && tar --zstd -C /usr -xf /tmp/ollama.tar.zst \
    && rm -f /tmp/ollama.tar.zst \
    && ollama --version

ENV OLLAMA_BASE_MODEL=${OLLAMA_BASE_MODEL}
RUN ollama serve & \
    OLLAMA_PID=$! && \
    sleep 3 && \
    curl -sf http://localhost:11434/ >/dev/null && \
    ollama pull "${OLLAMA_BASE_MODEL}" && \
    kill "${OLLAMA_PID}"

# This is just the name the pipeline will call at runtime -- the
# actual pfml-parser model tag doesn't exist yet. It's created by
# entrypoint.sh once your source is mounted in and readable.
ENV OLLAMA_MODEL=pfml-parser

# ============================================================
# 9. [REMOVED] Application code
#    No longer COPYed into the image -- bind-mount your working
#    directory to /app at `docker run` time instead:
#      docker run -it -v $(pwd):/app msbase
#    entrypoint.sh verifies the mount landed correctly before
#    proceeding, and bakes pfml-parser from the mounted source.
# ============================================================

# ============================================================
# 10. [REMOVED] Compile microsim_gp
#     Compiling against microsim_gp.c is now a manual step you run
#     yourself inside the container once your source is mounted,
#     e.g.:
#       mpicc -o /usr/local/bin/microsim_gp /app/microsim_gp/microsim_gp.c \
#           -I/app/microsim_gp/functions -I/app/microsim_gp/solverloop \
#           $(pkg-config --cflags hdf5-openmpi) $(pkg-config --libs hdf5-openmpi) \
#           -lm -lgsl -lgslcblas
#     The toolchain and libraries this needs (sections 3-5) are
#     already present in the image.
# ============================================================

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENV MICROSIM_DIR=/app/microsim_gp \
    MPI_NP=${MPI_NP} \
    MPI_GRID=${MPI_GRID}

WORKDIR /app
ENTRYPOINT ["entrypoint.sh"]
CMD ["bash"]

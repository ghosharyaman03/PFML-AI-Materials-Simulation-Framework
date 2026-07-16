FROM ubuntu:22.04

ARG OLLAMA_BASE_MODEL=gemma3:4b
ARG MPI_NP=4
ARG MPI_GRID="2 2"

ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PYTHONUNBUFFERED=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

SHELL ["/bin/bash", "-euo", "pipefail", "-c"]

# Reliable APT retries/timeouts
RUN echo 'Acquire::Retries "5";'              > /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::http::Timeout "30";'       >> /etc/apt/apt.conf.d/80-retries && \
    echo 'Acquire::https::Timeout "30";'      >> /etc/apt/apt.conf.d/80-retries

# Core tools; zstd is required for tar --zstd below.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        wget \
        zstd \
    && rm -rf /var/lib/apt/lists/*

# Compiler and build toolchain
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        gcc \
        g++ \
        gfortran \
        make \
        pkg-config \
        cmake \
    && rm -rf /var/lib/apt/lists/*

# MPI
RUN apt-get update && apt-get install -y --no-install-recommends \
        openmpi-bin \
        libopenmpi-dev \
    && rm -rf /var/lib/apt/lists/*

# Scientific libraries
RUN apt-get update && apt-get install -y --no-install-recommends \
        libhdf5-openmpi-dev \
        hdf5-tools \
        libgsl-dev \
        libgsl27 \
        libopenblas-dev \
        liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

# Python
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-dev \
        python-is-python3 \
    && rm -rf /var/lib/apt/lists/* \
    && pip3 install --no-cache-dir --upgrade pip

# Python dependencies
RUN pip3 install --no-cache-dir \
        numpy \
        pandas \
        scipy \
        pycalphad \
        h5py \
        matplotlib

RUN --mount=type=cache,target=/root/.cache/ollama-download,sharing=locked \
    curl -fL --retry 50 --retry-delay 5 --retry-all-errors --retry-connrefused \
        -C - --connect-timeout 30 \
        https://ollama.com/download/ollama-linux-amd64.tar.zst \
        -o /root/.cache/ollama-download/ollama.tar.zst \
    && tar --zstd -C /usr -xf /root/.cache/ollama-download/ollama.tar.zst \
    && ollama --version

# Download the base model at image-build time.
ENV OLLAMA_BASE_MODEL=${OLLAMA_BASE_MODEL}

RUN ollama serve & \
    OLLAMA_PID=$! && \
    sleep 3 && \
    curl -sf http://localhost:11434/ >/dev/null && \
    ollama pull "${OLLAMA_BASE_MODEL}" && \
    kill "${OLLAMA_PID}"

# entrypoint.sh builds this custom model from the mounted project source.
ENV OLLAMA_MODEL=pfml-parser

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENV MICROSIM_DIR=/app/microsim_gp \
    MPI_NP=${MPI_NP} \
    MPI_GRID=${MPI_GRID}

WORKDIR /app

ENTRYPOINT ["entrypoint.sh"]
CMD ["bash"]

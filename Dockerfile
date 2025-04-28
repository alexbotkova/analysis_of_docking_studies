FROM ubuntu:22.04

ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    wget \
    bzip2 \
    build-essential \
    libgl1-mesa-glx \
    libglu1-mesa \
    libxrender1 \
    libxext6 \
    libsm6 \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    bash /miniconda.sh -b -p $CONDA_DIR && \
    rm /miniconda.sh && \
    conda clean --all -y

RUN conda update -n base -c defaults conda -y

COPY requirements.txt .

RUN conda install -n base -c conda-forge python=3.10 pymol-open-source -y && \
    pip install -r requirements.txt

COPY . /app
WORKDIR /app

ENTRYPOINT ["/bin/bash"]
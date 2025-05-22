ARG BASE_IMAGE=bioconductor/bioconductor_docker:RELEASE_3_21
FROM ${BASE_IMAGE}

LABEL maintainer="Benjamin Haibe-Kains <benjamin.haibe.kains@utoronto.ca>"
LABEL description="Docker image for the PharmacoGx R/Bioconductor package"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
# Install pak
RUN Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) { install.packages("pak") }'

# Install package dependencies using pak - extract them from DESCRIPTION
RUN Rscript -e 'pak::pkg_install("bhklab/PharmacoGx", ask = FALSE, dependencies = TRUE, upgrade = FALSE)'

# Copy the local package files
COPY . /app

# Set working directory
WORKDIR /app

# Default command when the container starts
CMD ["R"]

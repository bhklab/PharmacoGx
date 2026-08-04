ARG BASE_IMAGE=bioconductor/bioconductor_docker:RELEASE_3_23
FROM ${BASE_IMAGE}

LABEL maintainer="Benjamin Haibe-Kains <benjamin.haibe.kains@utoronto.ca>"
LABEL description="Docker image for the PharmacoGx R/Bioconductor package"

# Package versions intentionally follow the selected Bioconductor base image.
# hadolint ignore=DL3008
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff-dev \
    libjpeg-dev \
    build-essential \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
# Install pak
RUN Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) { install.packages("pak") }'

# Install package dependencies separately so source-only changes can reuse this
# layer.
WORKDIR /app
COPY DESCRIPTION /app/DESCRIPTION
RUN Rscript -e 'pak::local_install_deps(".", ask = FALSE, dependencies = TRUE, upgrade = FALSE)'

# Copy the local package files
COPY . /app

# Install the checked-out source rather than the default GitHub branch.
RUN R CMD INSTALL --no-multiarch --with-keep.source .

# Default command when the container starts
CMD ["R"]

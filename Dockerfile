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
RUN Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) { \
    install.packages("pak", repos = "https://cloud.r-project.org") }'

# Copy only dependency information first to leverage Docker caching
COPY DESCRIPTION /tmp/

# Install package dependencies using pak - extract them from DESCRIPTION
RUN cd /tmp && \
    Rscript -e 'deps <- pak::pkg_deps_tree(".", dependencies = TRUE)$package; \
    pak::pkg_install(deps, ask = FALSE)'

# Copy the local package files
COPY . /app

# Set working directory
WORKDIR /app

# Build and install the package
RUN R CMD build . && \
    R CMD INSTALL *.tar.gz

# Default command when the container starts
CMD ["R"]

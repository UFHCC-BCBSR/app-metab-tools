# Use the official R Shiny image as the base
FROM rocker/shiny:latest

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    libfontconfig1-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages from CRAN
RUN install2.r --error \
    shiny \
    shinydashboard \
    DT \
    plotly \
    readxl \
    dplyr \
    tibble \
    shinyjs \
    BiocManager

# Install Bioconductor packages
RUN R -e "BiocManager::install(c('sva', 'limma'))"

# Expose the default Shiny port
EXPOSE 3838

# Run the Shiny app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server', host='0.0.0.0', port=3838)"]

FROM rocker/rstudio:4.1.0
LABEL Sam Widmayer <sjwidmay@gmail.com>

RUN  apt-get --allow-releaseinfo-change update \
    && apt-get install -y procps \
    ssh \
    bash \
    pkg-config \
    libglpk-dev \
    libjpeg62-dev \
    libz-dev \
    tk \
    libxml2 \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    xterm \
    x11-utils \
    libcairo2-dev \
    libblas-dev \
    libssh2-1-dev \
    libgit2-dev
    
RUN R -e "install.packages('purrr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('dplyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('tidyr', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('ggplot2', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('plotly', repos='http://cran.us.r-project.org')"
RUN R -e "install.packages('knitr', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('png', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('rmarkdown', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('pandoc', repos='http://cran.us.r-project.org', dependencies = TRUE)"
RUN R -e "install.packages('remotes', repos='http://cran.us.r-project.org')"
RUN R -e "remotes::install_github('rqtl/qtl2')"
RUN R -e "remotes::install_github('kbroman/broman')"
RUN R -e "remotes::install_github('kbroman/qtlcharts')"
RUN R -e "install.packages('cowplot')"

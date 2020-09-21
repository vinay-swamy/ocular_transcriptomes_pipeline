FROM rocker/verse:3.5.0
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    git 

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    zlib1g-dev

#STAR
RUN wget --quiet --no-check-certificate https://github.com/alexdobin/STAR/archive/2.6.1c.tar.gz && \
    tar -xzf 2.6.1c.tar.gz && \
    cd STAR-2.6.1c/source && \
    make STAR &&\
    echo DONE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    pkg-config \
    gawk \
    make \
    gcc \
    build-essential \
    libboost-all-dev \
    libz-dev libbz2-dev \
    autoconf \
    automake \
    libtool \
    unzip \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libncurses5-dev

#stringtie
RUN cd / && \
    wget --quiet --no-check-certificate https://github.com/gpertea/stringtie/archive/v2.0.3.tar.gz && \
    tar -xzf v2.0.3.tar.gz &&  \
    cd stringtie-2.0.3/ && \
    make release 

#htslib
RUN wget --quiet --no-check-certificate -O htslib.tar.gz https://github.com/samtools/htslib/archive/1.9.tar.gz && \
    tar -xzf htslib.tar.gz && \
    cd htslib-1.9 && \
    autoreconf && \
    ./configure && \
    make && \
    make install

#samtools
RUN wget --quiet --no-check-certificate https://github.com/samtools/samtools/archive/1.9.tar.gz && \
    tar -xzf 1.9.tar.gz && \
    cd samtools-1.9 && \ 
        autoreconf  && \
        ./configure && \
        make && \
        make install
#gffcompare
RUN wget --quiet --no-check-certificate https://github.com/gpertea/gffcompare/archive/v0.11.2.tar.gz && \
    tar -xzf v0.11.2.tar.gz && \
    git config --global http.sslVerify false && \
    git clone https://github.com/gpertea/gclib.git && \
    cd gffcompare-0.11.2/ && \
    make release

#bedops
RUN wget --quiet --no-check-certificate https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2 && \
    tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2 && \
    cp bin/*  /usr/local/bin

#salmon
RUN wget --quiet --no-check-certificate https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz && \
    tar -xzf salmon-1.1.0_linux_x86_64.tar.gz && \
    echo "export PATH=/$PWD/salmon-1.1.0_linux_x86_64/bin:${PATH}" >> /root/.bashrc

# conda 
RUN apt-get -qq -y install bzip2 curl \
    && mkdir  /conda \
    && curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /conda 
#VEP and python 
RUN /conda/bin/conda install -y python=3  \
    && /conda/bin/conda install -y -c bioconda ensembl-vep=100.4 \
    && /conda/bin/conda update conda \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && /conda/bin/conda clean --all --yes

#TransDecoder
RUN wget --quiet --no-check-certificate https://github.com/TransDecoder/TransDecoder/archive/v5.0.1.tar.gz && \
    tar -xzf v5.0.1.tar.gz && \
    echo "export PATH=/$PWD/TransDecoder-5.0.1:${PATH}" >> /root/.bashrc


#UCSC
RUN apt-get update \
    && apt-get install -y --no-install-recommends  rsync \
    && mkdir uscs/  \
    && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ ucsc/  \
    && echo "export PATH=/$PWD/ucsc:${PATH}" >> /root/.bashrc

#mosdepth
RUN wget --quiet --no-check-certificate -O /usr/local/bin/mosdepth  https://github.com/brentp/mosdepth/releases/download/v0.2.3/mosdepth && \
    chmod 777 /usr/local/bin/mosdepth

#hmmer 
RUN wget --quiet --no-check-certificate http://eddylab.org/software/hmmer/hmmer-3.2.1.tar.gz && \
    tar -xzf hmmer-3.2.1.tar.gz && \
    cd hmmer-3.2.1 && \
    ./configure && \
    make && \
    make check && \
    make install            
#bedtools
RUN wget --quiet --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.27.1/bedtools-2.27.1.tar.gz  && \
    tar -xzf bedtools-2.27.1.tar.gz && \
    cd bedtools2/ && \
    make && \
    make install && \
    bedtools --help

#R
COPY Rprofile.site /usr/local/lib/R/etc/Rprofile.site
RUN install2.r --deps TRUE \
    BiocManager \
    BiocVersion \
    renv


RUN echo ""
COPY renv.lock .
RUN apt-get update \
    && apt-get install -y --no-install-recommends libudunits2-dev && \
     Rscript -e "renv::init()"


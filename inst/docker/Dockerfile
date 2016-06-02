############################################################
# Dockerfile to build hisat2 container images
# Based on Ubuntu
############################################################

# Set the base image to Ubuntu
FROM bioconductor/devel

# File Author / Maintainer
MAINTAINER Sean Davis

# Configure default locale
RUN echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

# Update the repository sources list
RUN apt-get update  && apt-get -y install samtools=0.1.19-1 && apt-get install -y wget

# Install compiler and perl stuff
RUN apt-get install --yes build-essential git \
   build-essential \
   curl \
   python \ 
   python2.7-dev \
   python-pip \
   python-numpy \
   python-matplotlib 

# get NBCI Toolkit
WORKDIR /tmp
RUN wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/ngs/1.2.4/ngs-sdk.1.2.4-linux.tar.gz" && \
    tar zxfv ngs-sdk.1.2.4-linux.tar.gz

# Compile
#RUN make USE_SRA=1 NCBI_NGS_DIR=/tmp/ngs-sdk.1.2.3-linux NCBI_VDB_DIR=/tmp/ngs-sdk.1.2.3-linux
#RUN cp -p hisat2 hisat2-* /usr/local/bin

WORKDIR /tmp

RUN apt-get clean

# Set default working path
WORKDIR /root

# Default command to execute at startup of the container
CMD /bin/bash
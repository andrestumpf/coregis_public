FROM ubuntu:16.04
MAINTAINER Andre Stumpf andre.stumpf@unistra.fr
ENV HOME /root
RUN echo $HOME


# required system libraries
RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y install software-properties-common \
                       wget
                       git \
                       cmake \
                       make \
                       imagemagick \
                       libimage-exiftool-perl \
                       exiv2 \
                       libgeo-proj4-perl \
                       libx11-dev \
                       cmake-curses-gui \
                       freeglut3-dev \
                       libpcre3 \
                       libpcre3-dev \
                       libtiff5 \
                       libtiff5-dev \
                       libgdal-dev && \
    apt-get clean all

# Install the MicMac
WORKDIR $HOME
RUN git clone https://github.com/micmacIGN/micmac.git
WORKDIR $HOME/micmac
RUN mkdir build
WORKDIR $HOME/micmac/build
RUN cmake -DWITH_QT5=OFF -DBUILD_POISSON=OFF ../
RUN make install -j$(nproc)
ENV PATH "$PATH:$HOME/micmac/bin"

# Install the Orfeo toolbox
RUN add-apt-repository ppa:ubuntugis/ubuntugis-unstable \
    apt-get update \
    apt-get -y install otb-bin

# Install Anaconda
WORKDIR $HOME
RUN wget https://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh && \
    bash Anaconda3-5.0.1-Linux-x86_64.sh -b -p ~/anaconda && \
    rm Anaconda3-5.0.1-Linux-x86_64.sh
ENV PATH "$PATH:$HOME/anaconda/bin"
RUN conda config --add channels conda-forge && \
    conda update conda

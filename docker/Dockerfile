# Create an image with that looks like a naive user on the rhino servers.
#
# This Dockerfile will compile and install both Beast and
#
FROM r-base
MAINTAINER Chris Warth cwarth@fredhutch.org

ENV HOME /root
WORKDIR ${HOME}

ADD Dockerfile ${HOME}/Dockerfile

RUN apt-get update -qq
RUN apt-get dist-upgrade -y

# Install tools needed to compile Beast and Beagle.
# 
RUN apt-get install -y ant build-essential autoconf automake libtool subversion pkg-config

RUN apt-get install -y \
    python-dev \
    wget \
    openssh-client  \
    git  \
    pandoc  \
    lmodern  \
    openjdk-7-jdk  \
    openjdk-7-jre-headless  \
    cmake  \
    samtools \
    freetype* \
    pkg-config \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev

RUN wget -q https://bootstrap.pypa.io/get-pip.py
RUN python get-pip.py
RUN pip install --upgrade distribute

RUN pip install virtualenv

# # Make ssh dir
WORKDIR /root

# Make ssh dir
RUN mkdir /root/.ssh/

# Copy over private key, and set permissions
ADD id_rsa-github_readonly_user /root/.ssh/id_rsa
RUN chmod 700  /root/.ssh/id_rsa

# Create known_hosts
RUN touch /root/.ssh/known_hosts

# Add bitbuckets key
RUN ssh-keyscan -T 60 github.com >> /root/.ssh/known_hosts

# clone repo and remove private key
RUN git clone --recursive git@github.com:matsengrp/founder-inference-bakeoff.git
RUN rm -f /root/.ssh/id_rsa

ADD Dockerfile /root/Dockerfile

# Clone and compile Beagle from sources
RUN git clone --depth=1 https://github.com/beagle-dev/beagle-lib.git
WORKDIR ${HOME}/beagle-lib
RUN ./autogen.sh 
RUN ./configure --disable-sse --disable-march-native --prefix=/usr/local
RUN make install
#RUN rm -rf ${HOME}/beagle-lib

# Clone and install Beast from sources
WORKDIR ${HOME}
RUN git clone https://github.com/beast-dev/beast-mcmc.git
WORKDIR ${HOME}/beast-mcmc
RUN ant linux
RUN mv ${HOME}/beast-mcmc/release/Linux/BEASTv1.8.3pre ${HOME}/
#RUN rm -rf ${HOME}/beast-mcmc
# this leaves the compiled BEAST jar file in ${HOME}/BEASTv1.8.3pre/


# make beagle libraries and Beast script accessible
RUN echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >> ${HOME}/.bashrc
ENV LD_LIBRARY_PATH /usr/local/lib:${LD_LIBRARY_PATH}

# make beast scripts accessilbe on PATH
RUN echo 'export PATH=BEASTv1.8.3pre/bin:$PATH' >> ${HOME}/.bashrc

# Clean-up
RUN apt-get -y autoremove
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

CMD bash



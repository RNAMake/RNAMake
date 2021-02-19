FROM centos:8

RUN     dnf -y install scl-utils scl-utils-build git python3 'dnf-command(config-manager)' && \
        dnf config-manager --enable powertools && \
        dnf install -y  gcc-toolset-9
RUN     dnf install -y -q glibc-static


SHELL ["/usr/bin/scl", "enable", "gcc-toolset-9"]



# getting cmake
RUN curl -sSL https://github.com/Kitware/CMake/releases/download/v3.17.2/cmake-3.17.2-Linux-x86_64.tar.gz | tar -zx -C /usr/local --strip-components=1
# getting Ninja 
RUN git clone https://github.com/ninja-build/ninja.git && \
    cd ninja && \
    python3 configure.py --bootstrap && \
    cp ninja /usr/bin


RUN mkdir /RNAMake
ADD . /RNAMake
ENV RNAMAKE=/RNAMake
ENV X3DNA=$RNAMAKE/resources/x3dna/linux
RUN cd RNAMake/cmake/build/ && \
    python3 ../make_project.py -target linux && \
    cmake -G Ninja && \
    ninja

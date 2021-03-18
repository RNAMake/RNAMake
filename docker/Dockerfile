FROM centos:8 as build

# Get all build files
RUN     dnf -y install scl-utils scl-utils-build git python3 'dnf-command(config-manager)' && \
        dnf config-manager --enable powertools && \
        dnf install -y gcc-toolset-9 ninja-build cmake
RUN     dnf install -y -q glibc-static

SHELL ["/usr/bin/scl", "enable", "gcc-toolset-9"]

# Build RNAMake
RUN mkdir /RNAMake
ADD . /RNAMake
ENV RNAMAKE=/RNAMake
ENV X3DNA=$RNAMAKE/resources/x3dna/linux
RUN cd RNAMake/cmake/build/ && \
    python3 ../make_project.py -target linux && \
    cmake -G Ninja && \
    ninja && \
    chmod +x ./*

# Copy binaries and libraries to minimize image
FROM centos:8 as stage
COPY --from=build /RNAMake/cmake/build /RNAMake/bin
COPY --from=build /RNAMake/resources /RNAMake/resources
COPY --from=build /RNAMake/apps /RNAMake/apps
ENV RNAMAKE=/RNAMake
ENV X3DNA=$RNAMAKE/resources/x3dna/linux
ENV PATH=$PATH:/RNAMake/bin
WORKDIR /work

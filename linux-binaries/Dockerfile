FROM centos:6

RUN yum -y update && yum clean all

RUN     yum -y -q update && \
        yum -y -q install   git  gcc gcc-c++ \
        scl-utils centos-release-scl 
RUN     yum install -y devtoolset-7-gcc* \
        rh-python36
RUN     yum install -y -q glibc-static

SHELL ["/usr/bin/scl", "enable", "devtoolset-7","rh-python36"]

# getting cmake
RUN curl -sSL https://github.com/Kitware/CMake/releases/download/v3.17.2/cmake-3.17.2-Linux-x86_64.tar.gz | tar -zx -C /usr/local --strip-components=1
# getting Ninja 
RUN git clone https://github.com/ninja-build/ninja.git && \
    cd ninja && \
    ./configure.py --bootstrap && \
    cp ninja /usr/bin && \
    cd /

ADD build.sh /
RUN  source build.sh

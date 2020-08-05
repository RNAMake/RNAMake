FROM centos:6

RUN yum  install -y -q centos-release-scl && \
    yum install -y -q devtoolset-7 sclo-git25 rh-python36 glibc-static && \
    curl -sSL https://github.com/Kitware/CMake/releases/download/v3.17.2/cmake-3.17.2-Linux-x86_64.tar.gz | tar -zx -C /usr/local --strip-components=1

RUN yum  install -y -q zlib-static readline-static && mkdir -p /opt/src && \
    curl -sSL https://www.sqlite.org/2020/sqlite-autoconf-3310100.tar.gz | tar -zx -C /opt/src && \
    pushd /opt/src/sqlite-autoconf-3310100 && \
    scl enable devtoolset-7 rh-python36 sclo-git25 "./configure --prefix=/usr/local --enable-shared=no --enable-static=yes && make && make install" && \
    ldconfig

ADD  linux_build.sh /

RUN mkdir  -p    /result && \
    chmod +x   linux_build.sh && \
    scl enable devtoolset-7 rh-python36 sclo-git25 ./linux_build.sh

WORKDIR /
ENTRYPOINT cat /rnamake-binaries-linux.tar.gz

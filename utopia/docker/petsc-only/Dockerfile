FROM alpine:latest

ENV CC clang
ENV CXX clang++

RUN apk update
RUN apk upgrade
RUN apk add --no-cache \
            clang curl build-base gfortran perl cmake python linux-headers lapack-dev valgrind-dev gdb

#############################################################################

RUN cd /tmp/ && curl -L http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.9.tar.gz| tar xz

RUN cd /tmp/ && cd `ls | grep petsc` && \
    ./configure PETSC_ARCH=dorcker-opt \
                CC=$CC \
                CXX=$CXX \
                --prefix=/opt/petsc \
                --download-mpich \
                --download-superlu_dist \
                --download-metis \
                --download-parmetis \
                --download-mumps \
                --download-scalapack \
                --download-hypre

RUN cd /tmp/ && cd `ls | grep petsc` && \
    make && \
    make install

ENV PETSC_DIR /opt/petsc

#############################################################################

CMD ["/bin/sh"]

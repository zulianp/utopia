#FROM canesin/slepc:latest
#ENV PETSC_DIR /opt/petsc
#RUN apk add --no-cache curl build-base clang cmake git gfortran linux-headers 
#CMD ["/bin/sh"]

######################################################

FROM canesin/openblas:latest

ENV CC clang
ENV CXX clang++

# Compile PETSc with latest clang
RUN apk add --no-cache curl build-base gfortran clang perl cmake python linux-headers && \
    cd /tmp/ && \
    curl -L http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.9.tar.gz| tar xz && \
    cd `ls | grep petsc` && \
    ./configure PETSC_ARCH=dorcker-opt \
            CC=$CC \
            CXX=$CXX \
            --prefix=/opt/petsc \
            --known-64-bit-blas-indices=true \
            --with-blas-lapack-lib=/opt/OpenBLAS/lib/libopenblas.a \
            --download-mpich \
            --download-superlu_dist \
            --download-metis \
            --download-parmetis \
            --download-mumps \
            --download-scalapack \
            --download-hypre \
            --with-debugging=0 \
            COPTFLAGS='-O2 -march=ivybridge -mtune=ivybridge' \
            CXXOPTFLAGS='-O2 -march=ivybridge -mtune=ivybridge' \
            FOPTFLAGS='-O2 -march=ivybridge -mtune=ivybridge' && \
    make && \
    make install && \
    apk del curl build-base gfortran clang perl cmake python linux-headers && \
    rm -rf /usr/share/man/* /tmp/* /var/cache/apk/*

# Setup PETSc location
ENV PETSC_DIR /opt/petsc

CMD ["/bin/sh"]

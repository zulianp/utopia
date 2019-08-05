FROM alpine:latest

ENV CC clang
ENV CXX clang++

RUN apk update
RUN apk upgrade
RUN apk add --no-cache \
            clang curl build-base gfortran perl cmake python linux-headers lapack-dev valgrind-dev gdb git

#############################################################################

RUN cd /tmp/ && curl -L http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz | tar xz
RUN cd /tmp/ && cd `ls | grep mpich` && \
    ./configure --prefix=/opt/mpich

RUN cd /tmp/ && cd `ls | grep mpich` && make && make install

RUN git clone https://github.com/trilinos/Trilinos.git trilinos_cloned
RUN cd trilinos_cloned && mkdir build && cd build && \
    cmake .. \
      -DCMAKE_INSTALL_PREFIX=/opt/trilinos \
      -DTPL_ENABLE_MPI=ON \
      -DMPI_BASE_DIR=/opt/mpich \
      -DTrilinos_ENABLE_Tpetra=ON \
      -DTrilinos_ENABLE_TpetraCore=ON \
      -DTrilinos_ENABLE_Belos=ON \
      -DTrilinos_ENABLE_Amesos2=ON \
      -DTrilinos_ENABLE_Ifpack2=ON \
      -DTrilinos_ENABLE_MueLu=ON \
      -DTrilinos_ENABLE_NOX=ON  \
      -DAmesos2_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DIfpack2_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DKOKKOS_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DMueLu_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DRTOp_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DStratimikos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DTeuchos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DThyra_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DTpetraCore_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DTpetra_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DTrilinos_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON


RUN cd trilinos_cloned && cd build && make && make install

ENV TRILINOS_DIR /opt/trilinos

#############################################################################

CMD ["/bin/sh"]

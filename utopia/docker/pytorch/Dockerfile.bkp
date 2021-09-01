FROM frolvlad/alpine-miniconda3:python3.7

#################################################################################
RUN apk update
RUN apk upgrade
RUN apk add --no-cache \
    clang-dev curl build-base clang cmake git linux-headers m4 gdb bash \
    valgrind-dev \
    libexecinfo-dev \
    openssh  \
    lapack-dev blas-dev \
    libgfortran gfortran \
    perl perl-io-gzip perl-json \
    python3 \
    python3-dev \
    nano \
    swig \
    doxygen m4 autoconf automake libtool

RUN conda update -n base -c defaults conda
RUN conda install --yes pytorch torchvision torchaudio cpuonly -c pytorch

################################################################################
#                             MPICH installation
################################################################################
WORKDIR /tmp
RUN curl -L http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz | tar xz
RUN cd `ls | grep mpich` && ./configure --prefix=/opt/mpich
RUN cd `ls | grep mpich` && make && make install
ENV MPI_DIR /opt/mpich
ENV PATH $PATH:$MPI_DIR/bin

################################################################################
#                             PETSc installation
################################################################################
RUN git clone -b release https://gitlab.com/petsc/petsc.git petsc
WORKDIR /tmp/petsc
RUN ./configure PETSC_ARCH=dorcker-opt \
        --prefix=/opt/petsc \
        --with-mpi=1 \
        --with-mpi-dir=$MPI_DIR \
        --download-superlu_dist \
        --download-metis \
        --download-mumps \
        --download-scalapack \
        --download-hypre \
        --with-debugging=0
RUN make
RUN make install

ENV PETSC_DIR /opt/petsc

################################################################################
#                             Trilinos installation
################################################################################
WORKDIR /tmp
RUN git clone https://github.com/trilinos/Trilinos.git trilinos
WORKDIR /tmp/trilinos
RUN mkdir build && cd build && \
    cmake .. \
      -DCMAKE_INSTALL_PREFIX=/opt/trilinos \
      -DTPL_ENABLE_MPI=ON \
      -DMPI_BASE_DIR=$MPI_DIR \
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
      -DXpetra_ENABLE_EXPLICIT_INSTANTIATION=ON \
      -DCMAKE_CXX_STANDARD=14
RUN cd build && make && make install

ENV TRILINOS_DIR /opt/trilinos

################################################################################
#                             Enviroment variables
################################################################################
WORKDIR /
ENV MPI_DIR /opt/mpich
ENV LDFLAGS -L/opt/mpich/lib/:/opt/petsc/lib/:/opt/trilinos/lib/
ENV CPATH /opt/mpich/include/:/opt/petsc/include/:/opt/trilinos/include/
ENV LD_LIBRARY_PATH=/opt/trilinos/lib/:/opt/mpich/lib/:/opt/petsc/lib/
ENV PATH $PATH:$PETSC_DIR/bin
ENV PATH $PATH:$TRILINOS_DIR/bin
ENV CC clang
ENV CXX clang++


################################################################################
#                             Utopia installation
################################################################################

RUN git clone https://zulianp@bitbucket.org/zulianp/utopia.git
RUN cd utopia && git checkout vc_simd && git submodule update --init --recursive
RUN mkdir /utopia/utopia/build
RUN cd /utopia/utopia/build && cmake .. -DUTOPIA_SCRIPTING=ON \
    -DUTOPIA_ENABLE_TRILINOS=OFF \
    -DCMAKE_SHARED_LINKER_FLAGS:STRING="-shared" \
    -DBUILD_SHARED_LIBS=ON 

RUN cd /utopia/utopia/build && make -j

COPY ./hello_utopia_torch.py  /utopia/utopia/build/scripting/hello_utopia_torch.py
RUN cd /utopia/utopia/build/scripting && python hello_utopia_torch.py
CMD ["/bin/sh"]
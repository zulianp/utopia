FROM zulianp/mpich:latest

# ENV CC clang
# ENV CXX clang++

ENV MPI_DIR /opt/mpich
RUN echo "mpi binaries:"
RUN ls $MPI_DIR/bin | grep mpi
ENV PATH $PATH:$MPI_DIR/bin

RUN apk update
RUN apk upgrade
RUN apk add --no-cache \
            clang curl build-base gfortran perl cmake python linux-headers lapack-dev valgrind-dev gdb git gfortran m4


#############################################################################
############################## PETSC ########################################
#############################################################################

RUN git clone -b maint https://gitlab.com/petsc/petsc.git petsc_cloned

RUN cd petsc_cloned && \
    ./configure \
                --with-mpi-dir=$MPI_DIR \
                --prefix=/opt/petsc \
                --download-scalapack \
                --download-hypre \
                --with-cxx-dialect=C++11 \
                --download-superlu_dist \
                --download-superlu \
                --download-mumps \
                --download-hdf5 \
                --download-netcdf \
                --download-pnetcdf \
                --download-exodusii \
                --download-zlib \
                --download-triangle \
                --download-ctetgen \
                -with-debugging=0

RUN cd petsc_cloned && \
    make && \
    make install

ENV PETSC_DIR /opt/petsc
ENV LDFLAGS -L/opt/petsc/lib
ENV CPATH opt/petsc/include

#############################################################################
############################## SLEPC ########################################
#############################################################################

RUN git clone https://bitbucket.org/slepc/slepc slepc_cloned
RUN cd slepc_cloned && git checkout maint

# Dependencies for Arpack
RUN apk add --no-cache \
          autoconf automake libtool

ENV F77=/opt/mpich/bin/mpif77
RUN cd slepc_cloned && ./configure --prefix=/opt/slepc --download-arpack
RUN cd slepc_cloned && \
    make && \
    make install


ENV CPATH $CPATH:/opt/slepc/include

#############################################################################
############################## TRILINOS #####################################
#############################################################################

RUN git clone https://github.com/trilinos/Trilinos.git trilinos_cloned
RUN cd trilinos_cloned && mkdir build && cd build && \
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


RUN cd trilinos_cloned && cd build && make && make install


ENV TRILINOS_DIR /opt/trilinos
ENV LD_LIBRARY_PATH=/opt/petsc/lib/:/opt/slepc/lib

RUN echo 'export LD_LIBRARY_PATH='$LD_LIBRARY_PATH   > setup.sh
RUN echo 'export MOONOLITH_DIR='$MOONOLITH_DIR      >> setup.sh
RUN echo 'export TRILINOS_DIR='$TRILINOS_DIR        >> setup.sh
RUN env | grep PATH                                 >> setup.sh
RUN echo 'apk add bash'                             >> setup.sh
RUN echo 'apk add bash-completion'                  >> setup.sh
RUN echo 'apk add bash-doc'                         >> setup.sh
RUN echo 'apk add swig'                             >> setup.sh
RUN echo 'export PATH=$PATH:/opt/petsc/bin'         >> setup.sh

RUN cat setup.sh
RUN cat setup.sh                                    >> /root/.profile


#############################################################################
#############################################################################

CMD ["/bin/sh"]

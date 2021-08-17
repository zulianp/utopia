FROM centos:latest

# Install baseline
RUN yum -y upgrade
RUN yum -y install make
RUN yum -y install git cmake
RUN yum -y install clang-devel
RUN yum -y install python3-devel
RUN yum -y install swig
RUN yum -y install mpich-devel
RUN yum -y install m4 autoconf automake libtool
RUN yum -y install lapack blas
RUN yum -y install nano

# Install torch libs
RUN pip3 install torch torchvision torchaudio

ENV PATH=$PATH:/usr/lib64/mpich/bin/
ENV MPI_DIR=/usr/lib64/mpich
ENV LD_LIBRARY_PATH=/usr/lib64/

# # ################################################################################
# # #                             Utopia installation
# # ################################################################################

WORKDIR /
RUN git clone https://zulianp@bitbucket.org/zulianp/utopia.git
RUN cd utopia && git checkout monotone_mg_refactor && git submodule update --init --recursive
RUN mkdir /utopia/utopia/build
RUN cd /utopia/utopia/build && cmake .. \
    -DUTOPIA_ENABLE_SCRIPTING=ON \
    -DUTOPIA_DEPENDENCIES_DIR=/utopia_dependencies \
    -DBLAS_LIBRARIES=/usr/lib64/libblas.so.3 \
    -DLAPACK_LIBRARIES=/usr/lib64/liblapack.so.3

# Install petsc with utopia/cmake installer
RUN cd /utopia/utopia/build && make petsc
ENV PETSC_DIR=/utopia_dependencies/petsc

# Fix for bug in cmake/swig/numpy integration
RUN ln -s /usr/local/lib64/python3.6/site-packages/numpy/core/include/numpy/ /usr/include/numpy
RUN cd /utopia/utopia/build && cmake .. && make -j4
# RUN cd /utopia/utopia/build && make -j4 complete

COPY ./hello_utopia_torch.py  /utopia/utopia/build/scripting/hello_utopia_torch.py
RUN cd /utopia/utopia/build/scripting && python3 hello_utopia_torch.py
CMD ["/bin/sh"]
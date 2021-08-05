FROM centos:latest
RUN yum upgrade
RUN yum install make
RUN yum install git cmake
RUN yum install clang-devel
RUN yum install python3-devel
RUN git clone https://github.com/spack/spack.git
RUN . spack/share/spack/setup-env.sh
RUN spack install anaconda3

# RUN apk update
# RUN apk upgrade
# RUN apk add --no-cache \
#     clang-dev curl build-base clang cmake git linux-headers m4 gdb bash \
#     valgrind-dev \
#     libexecinfo-dev \
#     openssh  \
#     lapack-dev blas-dev \
#     libgfortran gfortran \
#     perl perl-io-gzip perl-json \
#     nano \
#     doxygen m4 autoconf automake libtool


# # COPY ./utopya_conda.yml ./utopya_conda.yml
# RUN conda update -n base -c defaults conda
# # RUN conda env update -n utopia_env --file ./utopya_conda.yml
# RUN conda init bash
# RUN source /root/.bashrc
# # RUN conda activate utopia_env
# # RUN conda install -y git
# # RUN conda install -y cmake
# # RUN conda install -y make
# # RUN conda install -y libgcc
# # RUN conda install -y openmpi openmpi-mpicc  openmpi-mpicxx  openmpi-mpifort

# WORKDIR /tmp
# RUN curl -L http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz | tar xz
# RUN cd `ls | grep mpich` && ./configure --prefix=/opt/mpich
# RUN cd `ls | grep mpich` && make && make install
# ENV MPI_DIR /opt/mpich
# ENV PATH $PATH:$MPI_DIR/bin


# # RUN find / -name "libstdc++.so"
# COPY ./mpi_test.cpp ./mpi_test.cpp
# RUN mpicxx mpi_test.cpp && mpiexec -n 1 ./a.out


# # Extra headers
# # RUN conda install -y -c conda-forge petsc
# # RUN conda install -y -c conda-forge trilinos
# # RUN conda install -y -c conda-forge doxygen
# # RUN conda install -y -c sarcasm clang-format
# # RUN conda install -y -c conda-forge libblas
# # RUN conda install -y -c conda-forge liblapack

# # ################################################################################
# # #                             Utopia installation
# # ################################################################################

# WORKDIR /
# RUN git clone https://zulianp@bitbucket.org/zulianp/utopia.git
# RUN cd utopia && git checkout development && git submodule update --init --recursive
# RUN mkdir /utopia/utopia/build
# RUN cd /utopia/utopia/build && cmake .. -DUTOPIA_DEPENDENCIES_DIR=/utopia_dependencies
# RUN cd /utopia/utopia/build && make petsc

# ENV PETSC_DIR=/utopia_dependencies/petsc
# # RUN cd /utopia/utopia/build && cmake .. && make -j4 complete


# RUN conda install --yes pytorch torchvision torchaudio cpuonly -c pytorch
# RUN conda install -y swig
# RUN conda install -y -c conda-forge python-devtools
# RUN conda install -y numpy-devel

# RUN mkdir /utopia/utopia/build_script
# RUN cd /utopia/utopia/build_script && cmake .. -DUTOPIA_SCRIPTING=ON

# # -DUTOPIA_SCRIPTING=ON \
# #     -DCMAKE_SHARED_LINKER_FLAGS:STRING="-shared" \
# #     -DBUILD_SHARED_LIBS=ON

# # RUN cd /utopia/utopia/build && make -j complete

# # COPY ./hello_utopia_torch.py  /utopia/utopia/build/scripting/hello_utopia_torch.py
# # RUN cd /utopia/utopia/build/scripting && python hello_utopia_torch.py
# CMD ["/bin/sh"]
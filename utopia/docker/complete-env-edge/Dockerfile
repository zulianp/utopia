FROM alpine:edge

# Add testing repository
RUN echo "http://dl-cdn.alpinelinux.org/alpine/edge/testing"     >> /etc/apk/repositories

RUN apk add --no-cache \
    clang-dev curl build-base clang cmake git linux-headers m4 gdb bash \
    valgrind-dev \
    libexecinfo-dev \
    openmpi openmpi-dev openssh  \
    lapack-dev blas-dev \
    libgfortran gfortran \
    perl perl-io-gzip perl-json \
    python2 \
    nano

#############################################################################
############################## LCOV ########################################
#############################################################################

RUN git clone https://github.com/linux-test-project/lcov.git /lcov_git
RUN cd /lcov_git && make install
RUN PATH=$PATH:/usr/local/bin/lcov

###################################################################################################################
###################################################################################################################

# Clone code
RUN git clone https://bitbucket.org/zulianp/utopia.git utopia_cloned
RUN cd /utopia_cloned && git checkout layout_primitive && git submodule update --init --recursive && git pull
RUN cd /utopia_cloned && git pull

###################################################################################################################
###################################################################################################################

# Compile petsc (dmplex dependencies do not compile, so we turn them off)
RUN cd /utopia_cloned/utopia && mkdir bin && cd bin && cmake .. -DUTOPIA_ENABLE_PETSC_DM_PLEX=OFF

#############################################################################
############################## PETSC ########################################
#############################################################################

RUN cd /utopia_cloned/utopia/bin && make petsc VERBOSE=1
ENV PETSC_DIR /utopia_cloned/utopia/external/petsc


#############################################################################
############################## TRILINOS #####################################
#############################################################################

RUN cd /utopia_cloned/utopia/bin && make trilinos VERBOSE=1
ENV TRILINOS_DIR /utopia_cloned/utopia/external/trilinos

#############################################################################
############################## SLEPC ########################################
#############################################################################

# RUN git clone https://bitbucket.org/slepc/slepc /slepc_cloned
# RUN cd /slepc_cloned && git checkout maint

# # Dependencies for Arpack
# RUN apk add --no-cache \
#           autoconf automake libtool

# ENV F77=mpif77
# RUN cd /slepc_cloned && ./configure --prefix=/utopia_cloned/utopia/external/slepc --download-arpack
# RUN cd /slepc_cloned && \
#     make && \
#     make install


# ENV SLEPC_DIR=/utopia_cloned/utopia/external/slepc
# ENV CPATH $CPATH:$SLEPC_DIR/include

###################################################################################################################
###################################################################################################################

ENV LD_LIBRARY_PATH=$PETSC_DIR/lib/:$TRILINOS_DIR/lib

RUN echo 'export LD_LIBRARY_PATH='$LD_LIBRARY_PATH   > setup.sh
RUN echo 'export PETSC_DIR='$PETSC_DIR              >> setup.sh
# RUN echo 'export SLEPC_DIR='$SLEPC_DIR              >> setup.sh
RUN echo 'export TRILINOS_DIR='$TRILINOS_DIR        >> setup.sh
RUN env | grep PATH                                 >> setup.sh
RUN echo 'apk add bash'                             >> setup.sh
RUN echo 'apk add bash-completion'                  >> setup.sh
RUN echo 'apk add bash-doc'                         >> setup.sh
RUN echo 'apk add swig'								>> setup.sh

RUN cat setup.sh
RUN cat setup.sh                                    >> /root/.profile

#############################################################################
#############################################################################

CMD ["/bin/sh"]

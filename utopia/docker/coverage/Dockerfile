FROM alpine:edge

# Add testing repository
RUN echo "http://dl-cdn.alpinelinux.org/alpine/edge/testing"     >> /etc/apk/repositories

RUN apk add --no-cache \
    clang-dev curl build-base clang cmake git linux-headers m4 gdb bash \
    valgrind-dev \
    openmpi openmpi-dev openssh  \
    lapack-dev blas-dev \
    netcdf-fortran libgfortran gfortran \
    perl perl-io-gzip perl-json \
    python2 \
    nano


RUN git clone https://github.com/linux-test-project/lcov.git /lcov_git
RUN cd /lcov_git && make install
RUN PATH=$PATH:/usr/local/bin/lcov

# Clone code
RUN git clone https://bitbucket.org/zulianp/utopia.git utopia_cloned
RUN cd /utopia_cloned && git checkout layout_primitive && git submodule update --init --recursive && git pull
RUN cd /utopia_cloned && git pull

# Compile petsc
RUN cd /utopia_cloned/utopia && mkdir bin && cd bin && cmake .. -DPETSC_CONFIG_ARGS="-hdf5-dir=/usr/"

# Compile
RUN cd /utopia_cloned/utopia && mkdir bin_cov
RUN cd /utopia_cloned/utopia/bin_cov && cmake -DENABLE_CODE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug ..
RUN cd /utopia_cloned/utopia/bin_cov && make complete

# Create coverage data
RUN cp /utopia_cloned/utopia/scripts/utils/coverage.sh /utopia_cloned/utopia/bin_cov/coverage.sh
RUN cd /utopia_cloned/utopia/bin_cov && source coverage.sh


RUN echo "For parallel runs use: mpirun --allow-run-as-root -n 2 ./utopia_test"

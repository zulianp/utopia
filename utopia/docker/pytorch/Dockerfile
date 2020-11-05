FROM frolvlad/alpine-miniconda3:python3.7

#################################################################################

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
    doxygen

RUN conda install --yes pytorch torchvision torchaudio cpuonly -c pytorch

#################################################################################
# TODO add missing installations (i.e., petsc and trilinos)
#################################################################################

RUN git clone https://zulianp@bitbucket.org/zulianp/utopia.git
RUN cd utopia && git checkout vc_simd && git submodule update --init --recursive
RUN mkdir /utopia/utopia/build
RUN cd /utopia/utopia/build && cmake .. -DUTOPIA_SCRIPTING=ON
RUN cd /utopia/utopia/build && make -j

COPY ./hello_utopia_torch.py  /utopia/utopia/build/scripting/hello_utopia_torch.py
RUN cd /utopia/utopia/build/scripting && python hello_utopia_torch.py
CMD ["/bin/sh"]
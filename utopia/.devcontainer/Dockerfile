# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.148.0/containers/cpp/.devcontainer/base.Dockerfile
# See tips for vscode: https://github.com/microsoft/vscode-tips-and-tricks

# [Choice] Debian / Ubuntu version: debian-10, debian-9, ubuntu-20.04, ubuntu-18.04
ARG VARIANT="buster"
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

# [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends git doxygen clang-format gfortran python3-dev swig cmake-curses-gui

RUN git clone https://github.com/spack/spack /home/vscode/spack && \
    . /home/vscode/spack/share/spack/setup-env.sh && export PATH=/home/vscode/spack/bin:$PATH && \
    spack install anaconda3 &&  \
    spack install openblas  &&  \
    spack install lapack    &&  \
    spack load anaconda3    &&  \
    conda install --yes  pytorch torchvision torchaudio cpuonly -c pytorch

# RUN /share/spack/setup-env.sh && spack install petsc +superlu-dist +mumps +metis +exodusii
# RUN . /home/vscode/spack/share/spack/setup-env.sh && spack install trilinos +tpetra +openmp +muelu +amesos2 +adios2 +ifpack2 +ml

RUN echo "#!/bin/bash" >                                                /home/vscode/init_container.sh
RUN echo "# Generate automatocally witin docker" >                      /home/vscode/init_container.sh
RUN echo "" >>                                                          /home/vscode/init_container.sh
RUN echo ". /home/vscode/spack/share/spack/setup-env.sh" >>             /home/vscode/init_container.sh
RUN echo "" >>                                                          /home/vscode/init_container.sh
RUN echo "export PATH=/home/vscode/spack/share/spack/bin:\$PATH"  >>    /home/vscode/init_container.sh
RUN echo "" >>                                                          /home/vscode/init_container.sh
RUN echo "spack load anaconda3; spack load openblas" >>                 /home/vscode/init_container.sh
RUN echo "" >>                                                          /home/vscode/init_container.sh
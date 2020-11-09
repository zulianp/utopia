# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.148.0/containers/cpp/.devcontainer/base.Dockerfile
# See tips for vscode: https://github.com/microsoft/vscode-tips-and-tricks

# [Choice] Debian / Ubuntu version: debian-10, debian-9, ubuntu-20.04, ubuntu-18.04
ARG VARIANT="bionic"
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

# # [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends git clang-format gfortran lcov python3-dev swig doxygen petsc-dev 
    #trilinos-all-dev (cmake is broken)


FROM centos:latest

# Install baseline
RUN yum upgrade
RUN yum install make
RUN yum install git cmake
RUN yum install clang-devel
RUN yum install python3-devel
RUN yum install swig
RUN yum install mpich-devel

# Install torch libs
RUN pip3 install torch torchvision torchaudio

ENV PATH=$PATH:/usr/lib64/mpich/bin/
ENV MPI_DIR=/usr/lib64/mpich

# # ################################################################################
# # #                             Utopia installation
# # ################################################################################

WORKDIR /
RUN git clone https://zulianp@bitbucket.org/zulianp/utopia.git
RUN cd utopia && git checkout monotone_mg_refactor && git submodule update --init --recursive
RUN mkdir /utopia/utopia/build
RUN cd /utopia/utopia/build && cmake .. -DUTOPIA_ENABLE_SCRIPTING=ON -DUTOPIA_DEPENDENCIES_DIR=/utopia_dependencies

# Install petsc with utopia/cmake installer
RUN cd /utopia/utopia/build && make petsc
ENV PETSC_DIR=/utopia_dependencies/petsc

RUN cd /utopia/utopia/build && cmake .. && make -j4 && make -j4 complete

COPY ./hello_utopia_torch.py  /utopia/utopia/build/scripting/hello_utopia_torch.py
RUN cd /utopia/utopia/build/scripting && python3 hello_utopia_torch.py
CMD ["/bin/sh"]
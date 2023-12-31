# utopia_fe/Dockerfile

# https://docs.docker.com/docker-hub/builds/automated-testing/
# https://medium.com/better-programming/how-to-test-your-automated-builds-on-docker-hub-e40879f35d1e
# https://github.com/NMichas/dockerbuild/blob/master/docker-compose.test.yml

FROM zulianp/utopia-fe-env:latest

ENV PETSC_DIR /opt/petsc
ENV MPI_DIR /opt/petsc
ENV LIBMESH_DIR /libmesh
ENV LDFLAGS -L/opt/petsc/lib
ENV CPATH opt/petsc/include:/opt/slepc/include

RUN ls $MPI_DIR/bin | grep mpi
RUN export PATH=$PATH:$MPI_DIR/bin

RUN apk add --no-cache curl build-base clang cmake git gfortran linux-headers m4 perl lapack-dev valgrind-dev gdb bash

# Docker hub stuff (comment if on PC)
COPY . /utopia_cloned
RUN ls /utopia_cloned


# RUN cd /utopia_cloned && git branch


# (Uncomment these ones if on PC)
# RUN git clone https://bitbucket.org/zulianp/utopia.git utopia_cloned
# # RUN cd /utopia_cloned && git checkout edsl_refactor && git submodule update --init --recursive && git pull
# RUN cd /utopia_cloned && git checkout development && git submodule update --init --recursive && git pull


##########################################################################
##########################################################################

# compile utopia

RUN cd /utopia_cloned/utopia && mkdir bin
RUN cd /utopia_cloned/utopia/bin && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/utopia -DCMAKE_BUILD_TYPE=Release

RUN cd /utopia_cloned/utopia/bin && make -j4
RUN cd /utopia_cloned/utopia/bin && make -j4 complete
RUN cd /utopia_cloned/utopia/bin && make install

##########################################################################

ENV LD_LIBRARY_PATH /opt/petsc/lib

# compile moonolith

RUN git clone https://zulianp@bitbucket.org/zulianp/par_moonolith.git moonolith_cloned
RUN cd /moonolith_cloned && mkdir bin
RUN cd /moonolith_cloned/bin && \
	cmake .. -DCMAKE_INSTALL_PREFIX=/par_moonolith -DCMAKE_BUILD_TYPE=Release \
    -DMPI_C_INCLUDE_PATH=$MPI_DIR/include -DMPI_C_LIBRARIES=$MPI_DIR/lib/libmpi.so  \
    -DMPI_CXX_LIBRARIES=$MPI_DIR/lib/libmpicxx.so -DMPI_CXX_INCLUDE_PATH=$MPI_DIR/include

RUN cd /moonolith_cloned/bin && make -j4 && make install


##########################################################################
##########################################################################

# compile utopia_fe

RUN cd /utopia_cloned/utopia_fe && mkdir bin
RUN cd /utopia_cloned/utopia_fe/bin && \
	cmake .. -DCMAKE_INSTALL_PREFIX=/utopia_fe -DUTOPIA_DIR=/utopia -DCMAKE_BUILD_TYPE=Release -DMOONOLITH_DIR=/par_moonolith -DUTOPIA_ENABLE_WARNINGS=OFF

RUN cd /utopia_cloned/utopia_fe/bin && make -j4
RUN cd /utopia_cloned/utopia_fe/bin && make -j4 complete
RUN cd /utopia_cloned/utopia_fe/bin && make install


##########################################################################
##########################################################################

ENV MOONOLITH_DIR=/par_moonolith
ENV UTOPIA_DIR /utopia
ENV UTOPIA_FE_DIR /utopia_fe
ENV LD_LIBRARY_PATH=/opt/petsc/lib/:/libmesh/lib/:/opt/slepc/lib

RUN echo 'export LD_LIBRARY_PATH='$LD_LIBRARY_PATH > setup.sh
RUN echo 'export MOONOLITH_DIR='$MOONOLITH_DIR    >> setup.sh
RUN echo 'export UTOPIA_FE_DIR='$UTOPIA_FE_DIR    >> setup.sh
RUN env | grep PATH 							  >> setup.sh
RUN echo 'apk add bash' 						  >> setup.sh
RUN echo 'apk add bash-completion'				  >> setup.sh
RUN echo 'apk add bash-doc'						  >> setup.sh
RUN echo 'export PATH=$PATH:/opt/petsc/bin'       >> setup.sh

RUN cat setup.sh
RUN cat setup.sh 								  >> /root/.profile

CMD ["/bin/sh"]

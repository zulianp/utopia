# utopia_fe/Dockerfile

FROM zulianp/utopia-trilinos-env:latest

ENV TRILINOS_DIR /opt/trilinos
ENV MPI_DIR /opt/mpich
ENV LDFLAGS -L/opt/mpich/lib 
#-L/opt/trilinos/lib
ENV CPATH /opt/mpich/include:/opt/trilinos/include
ENV LD_LIBRARY_PATH=/opt/trilinos/lib/:/opt/mpich/lib

RUN apk add --no-cache curl build-base clang cmake git gfortran linux-headers m4 perl lapack-dev valgrind-dev gdb

RUN git clone https://bitbucket.org/zulianp/utopia.git utopia_cloned
RUN cd /utopia_cloned && git checkout development && git submodule update --init --recursive && git pull

# compile utopia
RUN cd /utopia_cloned/utopia && mkdir bin 
RUN cd /utopia_cloned/utopia/bin && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/utopia -DCMAKE_BUILD_TYPE=Release  -DTRY_WITH_TRILINOS=ON


RUN cd /utopia_cloned/utopia/bin && make -j
RUN cd /utopia_cloned/utopia/bin && make install

RUN cd /utopia_cloned/utopia/bin && ./utopia_exec -verbose && ./utopia_exec -test none -bench

RUN echo 'export LD_LIBRARY_PATH='$LD_LIBRARY_PATH > setup.sh
RUN cat setup.sh

CMD ["/bin/sh"]

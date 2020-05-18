FROM zulianp/utopia-env:latest
ENV PETSC_DIR /opt/petsc
ENV LDFLAGS -L/opt/petsc/lib
ENV CPATH /opt/petsc/include:/opt/slepc/include

RUN apk add --no-cache curl build-base clang cmake git gfortran linux-headers m4 perl
RUN git clone https://github.com/libMesh/libmesh.git libmesh_git
RUN cd /libmesh_git && git submodule update --init
RUN cd /libmesh_git && ./configure  --prefix=/libmesh \
                                    --with-methods="opt dbg" \
                                    --disable-warnings \
                                    --enable-silent-rules \
                                    --enable-unique-id \
                                    --disable-warnings \
                                    --disable-maintainer-mode \
                                    --enable-petsc-hypre-required \
                                    --enable-metaphysicl-required \
                                    --with-mpi=/opt/petsc/

RUN cd /libmesh_git && make -j4
RUN cd /libmesh_git && make install
ENV LIBMESH_DIR /libmesh

CMD ["/bin/sh"]

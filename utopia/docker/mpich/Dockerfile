FROM alpine:latest

ENV CC clang
ENV CXX clang++

RUN apk update
RUN apk upgrade
RUN apk add --no-cache \
            clang curl build-base gfortran perl cmake python linux-headers lapack-dev valgrind-dev gdb git gfortran m4

#############################################################################

RUN cd /tmp/ && curl -L http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz | tar xz
RUN cd /tmp/ && cd `ls | grep mpich` && \
    ./configure --prefix=/opt/mpich

RUN cd /tmp/ && cd `ls | grep mpich` && make && make install

CMD ["/bin/sh"]
FROM canesin/slepc:latest

ENV PETSC_DIR /opt/petsc

# Enable when slepc wrapper is fixed to work without arpack
# ENV SLEPC_DIR /opt/slepc

RUN apk add --no-cache curl build-base clang cmake git gfortran linux-headers 

CMD ["/bin/sh"]

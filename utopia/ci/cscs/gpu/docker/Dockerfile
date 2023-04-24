ARG BASE_IMAGE
FROM $BASE_IMAGE

ARG NUM_PROCS

COPY . /utopia.src

RUN cd /opt/view/lib/pkgconfig && ln -s openblas.pc blas.pc

# scripting broken, hence we disable it.
# We are grepping in the output of cmake for the string 'Manually-specified variables were' to catch early error were variables are being renamed.
# It is an error to specify variables on the commandline that are unused.
RUN mkdir -p /utopia.src/build \
  && cd /utopia.src/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=/utopia.install \
    -DCMAKE_VERBOSE_MAKEFILE=ON \
    -DCMAKE_CXX_COMPILER=/opt/view/bin/nvcc_wrapper \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DBUILD_SHARED_LIBS=OFF \
    -DUTOPIA_ENABLE_SCRIPTING=OFF \
    -DUTOPIA_ENABLE_CXX14_FEATURES=ON \
    -DUTOPIA_ENABLE_PETSC=ON \
    -DUTOPIA_ENABLE_TRILINOS=ON \
    -DUTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE=ON \
    -DBLA_PREFER_PKGCONFIG=ON \
    ../utopia |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS complete \
  && make install

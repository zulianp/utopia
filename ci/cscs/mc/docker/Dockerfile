ARG BASE_IMAGE
FROM $BASE_IMAGE

ARG NUM_PROCS

COPY . /utopia.src

RUN cd /opt/view/lib/pkgconfig && ln -s openblas.pc blas.pc

ENV UTOPIA_DIR=/utopia.install

# scripting broken, hence we disable it.
# We are grepping in the output of cmake for the string 'Manually-specified variables were' to catch early error were variables are being renamed.
# It is an error to specify variables on the commandline that are unused.
RUN mkdir -p /utopia.src/build \
  && cd /utopia.src/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR \
    -DCMAKE_VERBOSE_MAKEFILE=OFF \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DCMAKE_CXX_STANDARD=17 \
    -DBUILD_SHARED_LIBS=OFF \
    -DUTOPIA_ENABLE_SCRIPTING=OFF \
    -DUTOPIA_ENABLE_CXX14_FEATURES=ON \
    -DUTOPIA_ENABLE_PETSC=ON \
    -DUTOPIA_ENABLE_TRILINOS=ON \
    -DUTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE=ON \
    -DBLA_PREFER_PKGCONFIG=ON \
    /utopia.src/utopia |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j$NUM_PROCS \
  && make install

RUN mkdir -p /utopia.src/utopia_fe/build \
  && cd /utopia.src/utopia_fe/build \
  && cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=$UTOPIA_DIR \
    -DCMAKE_PREFIX_PATH=$MARS_DIR/lib/cmake \
    -DCMAKE_VERBOSE_MAKEFILE=OFF \
    -DCMAKE_CXX_EXTENSIONS=OFF \
    -DCMAKE_CXX_STANDARD=17 \
    -DBUILD_SHARED_LIBS=OFF \
    -DUTOPIA_ENABLE_CXX14_FEATURES=ON \
    -DUTOPIA_ENABLE_LIBMESH=ON \
    -DUTOPIA_ENABLE_MOONOLITH=ON \
    -DUTOPIA_INSTALL_MOONOLITH=ON \
    -DUTOPIA_ENABLE_MARS=ON \
    /utopia.src/utopia_fe |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make par_moonolith \
  && cmake /utopia.src/utopia_fe  \
  && make -j$NUM_PROCS \
  && make install
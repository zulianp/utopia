FROM finkandreas/ubuntu-utopia-buildbase:20.04

COPY . /utopia.src

SHELL ["bash", "--rcfile", "/etc/profile", "-l", "-c"]
# scripting broken, hence we disable it.
# We are grepping in the output of cmake for the string 'Manually-specified variables were' to catch early error were variables are being renamed.
# It is an error to specify variables on the commandline that are unused.
RUN mkdir -p /utopia.src/build \
  && cd /utopia.src/build \
  && cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DBUILD_SHARED_LIBS=OFF  -DUTOPIA_ENABLE_SCRIPTING=OFF -DUTOPIA_REMOVE_TRILINOS_DEPRECATED_CODE=ON \
           -DUTOPIA_ENABLE_CXX14_FEATURES=ON -DCMAKE_CXX_STANDARD=17 -DCMAKE_CXX_EXTENSIONS=OFF \
           -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=/utopia.install ../utopia |& tee cmake_output.log \
  && ! grep 'Manually-specified variables were' cmake_output.log \
  && make -j4 complete \
  && make install


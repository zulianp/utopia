FROM finkandreas/ubuntu-utopia-buildbase:20.04

COPY . /utopia.src

SHELL ["bash", "--rcfile", "/etc/profile", "-l", "-c"]
RUN mkdir -p /utopia.src/build \
  && cd /utopia.src/build \
  && cmake -DCMAKE_VERBOSE_MAKEFILE=ON -DBUILD_SHARED_LIBS=OFF  -DUTOPIA_SCRIPTING=ON -DUTOPIA_ENABLE_CXX14_FEATURES=ON -DUTOPIA_ENABLE_PETSC=ON -DUTOPIA_ENABLE_TRILINOS=ON -DCMAKE_INSTALL_PREFIX=/utopia.install ../utopia \
  && make -j4 complete \
  && make install


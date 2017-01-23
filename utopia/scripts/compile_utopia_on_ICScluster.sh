#/bin/bash
dir=$PWD
utopia_dir="$(dirname "$dir")"
installation_suffix='_installation'

utopia_installation_dir=$utopia_dir$installation_suffix

cd ../bin 
rm -r *

cmake .. -DUTOPIA_STATIC_DEPENDENCIES_ONLY=ON -DUTOPIA_ARCHIVE_ONLY=ON -DPETSC_ARCH="" -DDOLFIN_SKIP_BUILD_TESTS=ON -DCMAKE_INSTALL_PREFIX=$utopia_installation_dir && make -j20 && make install

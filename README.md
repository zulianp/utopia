# Utopia #
Utopia is a C++ embedded domain specific language designed for parallel non-linear solution strategies and finite element analysis.

# Main contributors of utopia

- Dr. Patrick Zulian (Lead developer)
- Alena Kopanicakova (Linear and non-linear solvers, Fenics and MOOSE interoperability)
- Dr. Maria Giuseppina Chiara Nestola (Integration of parallel transfer for libmesh and MOOSE using Moonolith)

Developed at the Institute of Computational Science, USI, Lugano, Switzerland (https://www.ics.usi.ch).

# License
The software is realized with NO WARRANTY and it is licenzed under BSD 3-Clause license (https://opensource.org/licenses/BSD-3-Clause)

# Copyright
Copyright (c) 2015 Institute of Computational Science - USI Università della Svizzera Italiana, ETH-Z Eidgenössische Technische Hochschule Zürich

# Dependencies
- PETSc (https://www.mcs.anl.gov/petsc/), must be compiled with MUMPS enabled
- LibMesh for the experimental FE module (https://github.com/libMesh)

# Getting started

## Downloading utopia

Clone the repository and its submodules

- git clone --recurse-submodules https://bitbucket.org/zulianp/utopia.git

or for older git versions

- git clone https://bitbucket.org/zulianp/utopia.git
- cd utopia
- git submodule update --init --recursive

## Compiling utopia
Define the utopia path (you can also add it to your .bash_profile)
export UTOPIA\_DIR=<The absolute path of where you want to install utopia>

Go to the folder utopia/utopia:

- mkdir bin
- cd bin
- cmake .. -DCMAKE\_INSTALL\_PREFIX=$UTOPIA_DIR
- make
- make install.


## Compiling utopia_fe
After compiling utopia.

Define the utopia\_fe path (you can also add it to your .bash_profile)
export UTOPIA\_FE\_DIR=<The absolute path of where you want to install utopia\_fe>

You need a limesh installation. Define the libmesh install directory
export LIBMESH\_DIR=<The aboslute path of where you installed libmesh>

Go to the folder utopia/utopia\_fe:

- mkdir bin
- cd bin
- cmake .. -DUTOPIA\_DIR=$UTOPIA\_DIR -DLIBMESH_DIR=$LIBMESH_DIR -DCMAKE\_INSTALL\_PREFIX=$UTOPIA_FE_DIR -DMOONOLITH\_INSTALL\_PREFIX=$UTOPIA_FE_DIR
- make 
- make install


All the headers and binaries should be in the desired folder in the following form
- include
- lib
- bin
- config (here you can find useful configuration file for your cmake or make build system)

Setting MOONOLITH\_INSTALL\_PREFIX is optional. But if you want to delete the contect of the bin folder then it is required.

## Compiling your code with utopia and Makefile

If you are using utopia with 'make' you can use the utopia-config.makefile in the $UTOPIA\_DIR/config folder as shown
in the example in the file utopia/utopia/example\_usage\_of\_utopia/Makefile

If you are using utopia\_fe with 'make' you can use the utopia_fe_config.makefile in the $UTOPIA\_DIR/config folder as shown
in the example in the file utopia/utopia\_fe/example\_usage\_of\_utopia\_fe/Makefile

## Compiling your code with utopia_fe and CMake

If you are using utopia\_fe with 'cmake' you can use the utopia_fe_config.cmake in the $UTOPIA\_DIR/config folder as shown
in the example in the file utopia/utopia\_fe/example\_usage\_of\_utopia\_fe/cmake\_example/CMakeLists.txt

The FindUtopia.cmake and FindUtopiaFE.cmake can be found in utopia/utopia\_fe/cmake and they can be copied to your project when needed.

## Compiling utopia on Windows

The suggested environment for compiling utopia on Windows is Cygwin.

By using the Cygwin setup utility, install the following packages and their dependencies:
- `gcc-core`
- `gcc-fortran`
- `gcc-g++`
- `gdb`
- `make`
- `cmake`
- `openmpi`
- `libopenmpi-devel`
- `libopenblas`
- `liblapack-devel`

If you want to use PETSc, we suggest to disable its shared libraries to avoid having to tinker with the PATH environment variable. Here is an example configuration:
```
./configure --download-mumps --download-scalapack --with-shared-libraries=0
```
You can also specify a custom install directory for PETSc by using `--prefix`, as normal.

Follow the steps above (Compiling utopia) to compile utopia itself.

## Extras
- -DUTOPIA_ARCHIVE_ONLY=ON Allows to compile only the archive file utopia.a (or .lib for windows)
- -DUTOPIA_STATIC_DEPENDENCIES_ONLY=ON Allows to restrict the linking to static libraries
- -DMPI_CXX_COMPILER=<the desired compiler>. Allows to set the mpi compiler for C++.


The API documentation of Utopia can be generated through Doxygen by using the command *make docs* from the bin folder after calling cmake ... The API documentation is generated in the *utopia/doc/api* folder in both HTML (see html/index.html) and LateX (see refman.tex) formats.


## GPU with Trilinos TPetra
Following flags need to be set before lauching your code with Trilinos on GPU
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export CUDA_LAUNCH_BLOCKING=1

## Multinode GPU jobs
Multinode GPU can be launched by exporting following flag
export CRAY_CUDA_MPS=1


## CMake users

a FindUtopia.cmake  and a FindUtopiaFE.cmake are available in the utopia/utopia_fe/cmake folder. Define UTOPIA_DIR in your shell environment and the use the cmake find_package.

## Contact

Join us on slack:
https://join.slack.com/t/ics-utopia/signup


## Citing Utopia
If you use Utopia for your research you can use the following bibtex entry to cite us

```bibtex
@misc{utopiagit,
	author = {Patrick Zulian and Alena Kopani{\v c}{\'a}kov{\'a} and Maria Chiara Giuseppina Nestola and Andreas Fink and Nur Fadel and Victor Magri and Teseo Schneider and Eric Botter},
	title = {{U}topia: {A} {C}++ embedded domain specific language for scientific computing. {G}it repository},
	url = {https://bitbucket.org/zulianp/utopia},
	howpublished = {https://bitbucket.org/zulianp/utopia},
	year = {2016}
}
```

## Docker containers
A docker container based on Alpine-Linux can be be found at https://hub.docker.com/r/utopiadev/utopia and downloaded with `docker pull utopiadev/utopia`. For the moment only `utopia-petsc` is supported for this container. 
You can use `docker image ls` to find the image and run it with `docker run -it <image_hash>`

Download docker at https://www.docker.com/products/docker-desktop


# More details coming soon!


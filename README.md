# Utopia #
Utopia is a C++ embedded domain specific language designed for parallel non-linear solution strategies and finite element analysis.

# Contributors of utopia

- Patrick Zulian (Lead developer)
- Alena Kopanicakova (Linear and non-linear solvers, Fenics and MOOSE interoperability)
- Dr. Maria Giuseppina Chiara Nestola (Integration of parallel transfer for libmesh and MOOSE using Moonolith) 
- Teseo Schneider (Visualization tool)

# Dependencies
- PETSc (https://www.mcs.anl.gov/petsc/)
- LibMesh for the experimental FE module (https://github.com/libMesh)

# Getting started

## Compiling utopia
Define the utopia path (you can also add it to your .bash_profile)
export UTOPIA\_DIR=<The aboslute path of where you want to install utopia>

Go to the folder utopia/utopia:

- mkdir bin
- cd bin
- cmake .. -DCMAKE\_INSTALL\_PREFIX=$UTOPIA_DIR
- make
- make install.


## Compiling utopia_fe
After compiling utopia

You need a limesh installation. Define the libmesh install directory
export LIBMESH\_DIR=<The aboslute path of where you installed libmesh>

Go to the folder utopia/utopia\_fe:

- mkdir bin
- cd bin
- cmake -DUTOPIA\_DIR=$UTOPIA\_DIR -DLIBMESH_DIR=$LIBMESH_DIR -DCMAKE\_INSTALL\_PREFIX=$UTOPIA_DIR -DMOONOLITH\_INSTALL\_PREFIX=$UTOPIA_DIR
- make 
- make install


All the headers and binaries should be in the desired folder in the following form
- include
- lib
- bin
- config (here you can find useful configuration file for your cmake or make build system)

Setting MOONOLITH\_INSTALL\_PREFIX is optional. But if you want to delete the contect of the bin folder then it is required.

## Compiling your code with utopia

If you are using utopia with 'make' you can use the utopia_config.makefile in the $UTOPIA\_DIR/config folder as shown
in the example in the file utopia/utopia/example\_usage\_of\_utopia/Makefile

If you are using utopia\_fe with 'make' you can use the utopia_fe_config.makefile in the $UTOPIA\_DIR/config folder as shown
in the example in the file utopia/utopia_fe/example\_usage\_of\_utopia\_fe/Makefile


- -DUTOPIA_ARCHIVE_ONLY=ON Allows to compile only the archive file utopia.a (or .lib for windows)
- -DUTOPIA_STATIC_DEPENDENCIES_ONLY=ON Allows to restrict the linking to static libraries
- -DMPI_CXX_COMPILER=<the desired compiler>. Allows to set the mpi compiler for C++.


The API documentation of Utopia can be generated through Doxygen by using the command *make docs* from the bin folder after calling cmake ... The API documentation is generated in the *utopia/doc/api* folder in both HTML (see html/index.html) and LateX (see refman.tex) formats.



# More details coming soon!
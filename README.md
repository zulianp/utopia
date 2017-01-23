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
To compile utopia, go to the folder utopia

- mkdir bin
- cd bin
- cmake .. -DCMAKE\_INSTALL\_PREFIX=<The aboslute path of where you want to install utopia>.
- make
- make install.

All the headers and binaries should be in the desired folder in the following form

- include
- lib
- bin
- config (here you can find useful configuration file for your cmake or make build system)


- -DUTOPIA_ARCHIVE_ONLY=ON Allows to compile only the archive file utopia.a (or .lib for windows)
- -DUTOPIA_STATIC_DEPENDENCIES_ONLY=ON Allows to restrict the linking to static libraries
- -DMPI_CXX_COMPILER=<the desired compiler>. Allows to set the mpi compiler for C++.


The API documentation of Utopia can be generated through Doxygen by using the command *make docs* from the bin folder after calling cmake ... The API documentation is generated in the *utopia/doc/api* folder in both HTML (see html/index.html) and LateX (see refman.tex) formats.



# More details coming soon!
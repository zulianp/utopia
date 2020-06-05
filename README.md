[![CodeFactor](https://www.codefactor.io/repository/bitbucket/zulianp/utopia/badge)](https://www.codefactor.io/repository/bitbucket/zulianp/utopia)

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

([More in-depth guide here](https://bitbucket.org/zulianp/utopia/wiki/Getting%20started))

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


## Contact

Join us on slack:
https://join.slack.com/t/ics-utopia/signup


## Citing Utopia
If you use Utopia for your research please use the following bibtex entry (or equivalent) to cite us

```bibtex
@misc{utopiagit,
	author = {Patrick Zulian and Alena Kopani{\v c}{\'a}kov{\'a} and Maria Chiara Giuseppina Nestola and Andreas Fink and Nur Fadel and Alessandro Rigazzi and Victor Magri and Teseo Schneider and Eric Botter and Jan Mankau and Rolf Krause},
	title = {{U}topia: {A} {C}++ embedded domain specific language for scientific computing. {G}it repository},
	url = {https://bitbucket.org/zulianp/utopia},
	howpublished = {https://bitbucket.org/zulianp/utopia},
	year = {2016}
}
```

Several components are the outcome of specific authors' research work.  For citing individual components add `--citations` to your run. For instance, 
```
./my_app --citations
```

## Docker containers
A docker container based on Alpine-Linux can be be found at https://hub.docker.com/r/utopiadev/utopia and downloaded with `docker pull utopiadev/utopia`. For the moment only `utopia-petsc` is supported for this container. 
You can use `docker image ls` to find the image and run it with `docker run -v <host_file_directory>:<image_file_directory> -it <image_hash>`, for instance  `docker run -v ~/Desktop/my_mesh_files:/my_mesh_files -it <image_hash>`

Download docker at https://www.docker.com/products/docker-desktop


# More details coming soon!


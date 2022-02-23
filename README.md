[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![CodeFactor](https://www.codefactor.io/repository/bitbucket/zulianp/utopia/badge)](https://www.codefactor.io/repository/bitbucket/zulianp/utopia)


# Utopia #
Utopia is a C++ embedded domain specific language designed for parallel non-linear solution strategies and finite element analysis.

# Main contributors of utopia

- Dr. Patrick Zulian (Lead developer, ICS)
- Dr. Alena Kopanicakova (Linear and nonlinear solvers, ICS)
- Dr. Maria Giuseppina Chiara Nestola (MOOSE integrations, ICS)
- Dr. Nur Aiman Fadel (CSCS)
- Andreas Fink (CSCS)

Developed at the Institute of Computational Science, USI, Lugano, Switzerland (https://www.ics.usi.ch).

# License
The software is realized with NO WARRANTY and it is licenzed under [BSD 3-Clause license](https://opensource.org/licenses/BSD-3-Clause)

# Copyright
Copyright (c) 2015 Institute of Computational Science - USI Università della Svizzera Italiana, ETH-Z Eidgenössische Technische Hochschule Zürich

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

For compiling the very experimental FE library follow the instructions [HERE](https://bitbucket.org/zulianp/utopia/wiki/Utopia%20FE).
In the future the Utopia FE library will be removed from this repository and will have its own repo (after a large scale refactoring).

### MacOS Big Sur
For avoiding problems with cmake define the following env variable `export SYSTEM_VERSION_COMPAT=1`.

## Contact

- Join us on [slack](https://join.slack.com/t/ics-utopia/signup)
- Join us or contact us through our [mailing list](https://lists.usi.ch/mailman3/postorius/lists/utopia-users.lists.usi.ch)


## Citing Utopia
If you use Utopia for your research please use the following bibtex entry (or equivalent) to cite us

```bibtex
@article{utopia2021,
    author = {Zulian, Patrick and Kopani{\v{c}}{\'{a}}kov{\'{a}}, Alena and Nestola, Maria G C and Fadel, Nur and Fink, Andreas and VandeVondele, Joost and Krause, Rolf},
    title = {Large scale simulation of pressure induced phase‐field fracture propagation using {U}topia},
    journal = {CCF Transactions on High Performance Computing},
    year = {2021},
    month = {06},
    abstract = {Non-linear phase field models are increasingly used for the simulation of fracture propagation problems. 
	The numerical simulation of fracture networks of realistic size requires the efficient parallel solution of large coupled non-linear systems. 
	Although in principle efficient iterative multi-level methods for these types of problems are available, they are not widely used in practice due to the complexity of their parallel implementation. 
	Here, we present Utopia, which is an open-source C++ library for parallel non-linear multilevel solution strategies. 
	Utopia provides the advantages of high-level programming interfaces while at the same time a framework to access low-level data-structures without breaking code encapsulation. 
	Complex numerical procedures can be expressed with few lines of code, and evaluated by different implementations, libraries, or computing hardware. 
	In this paper, we investigate the parallel performance of our implementation of the recursive multilevel trust-region (RMTR) method based on the Utopia library. 
	RMTR is a globally convergent multilevel solution strategy designed to solve non-convex constrained minimization problems. 
	In particular, we solve pressure-induced phase-field fracture propagation in large and complex fracture networks. 
	Solving such problems is deemed challenging even for a few fractures, however, here we are considering networks of realistic size with up to 1000 fractures.},
    doi = {10.1007/s42514-021-00069-6},
    url = {https://doi.org/10.1007/s42514-021-00069-6},
    eprint = {https://doi.org/10.1007/s42514-021-00069-6},
}


@misc{utopiagit,
	author = {Patrick Zulian and Alena Kopani{\v c}{\'a}kov{\'a} and Maria Chiara Giuseppina Nestola and Andreas Fink and Nur Fadel and Alessandro Rigazzi and Victor Magri and Teseo Schneider and Eric Botter and Jan Mankau and Rolf Krause},
	title = {{U}topia: A performance portable {C}++ library for parallel linear and nonlinear algebra. {G}it repository},
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
You can use a pre-installed version of utopia in a Docker container. See [HERE](https://bitbucket.org/zulianp/utopia/wiki/Docker%20containers) for more details.


# More details coming soon!


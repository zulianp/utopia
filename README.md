[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![CodeFactor](https://www.codefactor.io/repository/bitbucket/zulianp/utopia/badge)](https://www.codefactor.io/repository/bitbucket/zulianp/utopia)


# Utopia #
Utopia is a C++ embedded domain specific language designed for parallel non-linear solution strategies and finite element analysis.

# Main contributors of utopia

- Dr. Patrick Zulian (Lead developer, Euler institute, USI)
- Dr. Alena Kopanicakova (Linear and nonlinear solvers, Euler institute, USI)
- Dr. Maria Giuseppina Chiara Nestola (MOOSE integrations, Euler institute, USI)
- Andreas Fink (CSCS)
- Dylan Ramelli (USI)
- Dr. Nur Aiman Fadel (CSCS)

Developed at the Euler institute, USI, Lugano, Switzerland (https://www.euler.usi.ch/).

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

## Contact

- Join us on [slack](https://join.slack.com/t/ics-utopia/signup)
- Join us or contact us through our [mailing list](https://lists.usi.ch/mailman3/postorius/lists/utopia-users.lists.usi.ch)


## Using Utopia for your research? Here is how you can cite it
If you use Utopia for your research please use the following bibtex entry to cite us

```bibtex
@article{utopia2021,
    author = {Zulian, Patrick and Kopani{\v{c}}{\'{a}}kov{\'{a}}, Alena and Nestola, Maria G C and Fadel, Nur and Fink, Andreas and VandeVondele, Joost and Krause, Rolf},
    title = {Large scale simulation of pressure induced phase‐field fracture propagation using {U}topia},
    journal = {CCF Transactions on High Performance Computing},
    year = {2021},
    month = {06},
    doi = {10.1007/s42514-021-00069-6},
    url = {https://doi.org/10.1007/s42514-021-00069-6},
    eprint = {https://doi.org/10.1007/s42514-021-00069-6},
}
```

Use the following bibtex to refer to the repository
```bibtex
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



## Contributions 

### Utopia Open-source repository (BSD 3-clause license)

1. https://bitbucket.org/zulianp/utopia
2. https://bitbucket.org/zulianp/par_moonolith (frontend and several integrations are in Utopia)


###  Simulations realized with Utopia


1. [Large-scale simulation of pressure-induced phase-field fracture propagation using Utopia](https://doi.org/10.1007/s42514-021-00069-6)
2. [3D non-conforming mesh model for flow in fractured porous media using Lagrange multipliers](https://doi.org/10.1016/j.cageo.2019.06.014)
3. [Comparison and application of non-conforming mesh models for flow in fractured porous media using dual Lagrange multipliers](https://doi.org/10.1016/j.jcp.2021.110773)

###  Contributions to benchmarks


1. [Verification benchmarks for single-phase flow in three-dimensional fractured porous media](https://doi.org/10.1016/j.advwatres.2020.103759)

###  Simulations using Utopia as a dependency

1. [Fully coupled dynamic simulations of bioprosthetic aortic valves based on an embedded strategy for fluid-structure interaction with contact](https://doi.org/10.1093/europace/euaa398)
2. [An immersed boundary method for fluid-structure interaction based on variational transfer](https://doi.org/10.1016/j.jcp.2019.108884)
3. [Space-time multilevel Monte Carlo methods and their application to cardiac electrophysiology](https://doi.org/10.1016/j.jcp.2021.110164)


# More details coming soon!


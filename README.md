[![CodeFactor](https://www.codefactor.io/repository/bitbucket/zulianp/utopia/badge)](https://www.codefactor.io/repository/bitbucket/zulianp/utopia) [![codecov](https://codecov.io/bb/zulianp/utopia/branch/master/graph/badge.svg)](https://codecov.io/bb/zulianp/utopia)

# Utopia #
Utopia is a C++ embedded domain specific language designed for parallel non-linear solution strategies and finite element analysis.

# Main contributors of utopia

- Dr. Patrick Zulian (Lead developer, ICS)
- Alena Kopanicakova (Linear and nonlinear solvers, ICS)
- Dr. Maria Giuseppina Chiara Nestola (MOOSE integrations, ICS)
- Nur Aiman Fadel (CSCS)
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

For compiling the very experimental FE library follow the instructions [HERE](https://bitbucket.org/zulianp/utopia/wiki/Utopia%20FE)

## Contact

- Join us on [slack](https://join.slack.com/t/ics-utopia/signup)
- Join us or contact us through our [mailing list](https://mailman2.ti-edu.ch/mailman/listinfo/utopia-users)


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
You can use a pre-installed version of utopia in a Docker container. See [HERE](https://bitbucket.org/zulianp/utopia/wiki/Docker%20containers) for more details.


# More details coming soon!


## About
utopia_compile.sh is an easy to use script for quickly setting up utopia. There are four different types of builds supported:
- basic: Build utopia with basic functionalities and with blas as backend.
- all: Build utopia with complete functionalities leveraginf PETSC and Trilinos as backend. This case assumes that you have these two dependencies already installed somewhere on your machine.
- fluya: Builds utopia in fluya mode.
- local: Builds petsc and trilinos along with utopia in a local folder inside your main utopia folder, i.e "utopia/external/petsc" and "utopia/external/Trilinos" respectively.


## Instructions on how to use utopia_compile.sh
While in "utopia/utopia/"
```
./scripts/utopia_compile.sh -b <build_type> -j <n_jobs> -p <install_prefix>
./scripts/utopia_compile.sh -h 
```
Where:

- build_type: basic, all, fluya, local.
- n_jobs: number of jobs you want to run.
- install_prefix: path to where you want to install utopia.
- -h: Will print out info on how to use the script as well.

## What the script actually does.
- The script will create a build_<build_type> folder where to compile utopia. If local selected it will also build and install petsc and trilinos.
- Will run utopia_bench and utopia_test.
- Will install utopia at <install_prefix> location.
- Finally will test the installation by creating a folder inside build_<build_type> and try to find utopia and use it to compile some test executables.

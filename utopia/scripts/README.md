# Build configurations

Utopia supports build configurations for the different applications

- FLUYA 	`-DUTOPIA_ENABLE_FLUYA_MODE=ON`
- AVFLOW	`-DUTOPIA_ENABLE_AVFLOW_MODE=ON`
- FRANETG 	`-DUTOPIA_ENABLE_FRANETG_MODE=ON`

Utopia provides the environment scripts for the supported supercomputers. Source them to load the relevant modules.

- Eiger: `source utopia/scripts/setup_fluya_env_eiger.sh`
- Daint: `source utopia/scripts/setup_fluya_env_daint.sh`

For compiling the main dependencies

## FLUYA

Some of the main dependencies can be compiled with utopia directly. Trilinos needs to be compiled separately (see `utopia/scripts/installer/configure_trilinos.sh`).
After running the environement script, go to `utopia` folder

```bash
mkdir -p build_fluya  && cd build_fluya
cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_INSTALL_PETSC=ON -DUTOPIA_PETSC_ENABLE_SUPERLU=ON -DUTOPIA_INSTALL_YAML_CPP=ON -DTrilinos_DIR=$TRILINOS_DIR/lib64/cmake/Trilinos -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fluya
make yaml-cpp
make petsc
cmake ..
make -j8 && make install
```

UtopiaFE is compiled with the following procedure. Go to `utopia_fe` folder

```bash
mkdir -p build_fluya &&
cd build_fluya &&
cmake .. -DUTOPIA_ENABLE_FLUYA_MODE=ON -DUTOPIA_INSTALL_MOONOLITH=ON -DUtopia_DIR=$INSTALL_DIR/utopia_fluya -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR/utopia_fe_fluya &&
make -j8 && make install
```

## AVFLOW

Similar but does not need Trilinos. LIBMESH can be installed directly too using `-DUTOPIA_INSTALL_LIBMESH` and `make libmesh`.


# Build scripts (experimental)
utopia_compile.sh is an easy to use script for quickly setting up utopia. There are four different types of builds supported:

- basic: Build utopia with basic functionalities and with blas as backend.
- all: Build utopia with complete functionalities leveraginf PETSC and Trilinos as backend. This case assumes that you have these two dependencies already installed somewhere on your machine.
- fluya: Builds utopia in fluya mode.
- local: Builds petsc and trilinos along with utopia in a local folder inside your main utopia folder, i.e "utopia/external/petsc" and "utopia/external/Trilinos" respectively.


## Instructions on how to use utopia_compile.sh
While in "utopia/utopia/"

```bash
./scripts/utopia_compile.sh -b <build_type> -j <n_jobs> -p <install_prefix>
./scripts/utopia_compile.sh -h
```
Where:

- build_type: basic, all, fluya, local.
- n_jobs: number of jobs you want to run.
- install_prefix: path to where you want to install utopia. Default installation folder is `/usr/local/`
- -h: Will print out info on how to use the script as well.

## What the script actually does.
- The script will create a build_<build_type> folder where to compile utopia. If local selected it will also build and install petsc and trilinos.
- Will run a cmake configuration based on the type of build selected. Before proceeding to the next step it will prompt the user to ask if the configuration was correct.
- Will run utopia_bench and utopia_test.
- Will install utopia at <install_prefix> location.
- Finally will test the installation by creating a folder inside build_<build_type> and try to find utopia and use it to compile some test executables.

## Extra
You can also run the different build scripts individually. For example:

```bash
./scripts/make_all.sh <n_jobs> <install_prefix>
```

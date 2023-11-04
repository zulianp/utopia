#!/usr/bin/env bash
set -e

# To crash if there is an error.
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

function display_help() {
    printf "This is the general utopia build script:\n ARGUMENTS:\n\t-h: Displays this information.\n\t-b <build_type>: Run the preferred build type.\n\t\t<build_type>: \n\t\tbasic: Compile utopia with only blas, lapack, umfpack\n\t\tall: Compile utopia with petsc, trilinos, blas, lapack, umfpack and other options.\n\t\tfluya: Compile utopia in fluya mode.\n\t\tlocal: Compile utopia and download/install petsc and trilinos locally in ../external/.\n\t-j <n_jobs>: Number of jobs to run in make.\n\t-p <install_prefix>: Path for utopia installation.\n"
}

function basic_build() {
    echo "Running basic build"
    ${SCRIPT_DIR}/make_basic.sh $1 $2
}


function all_build() {
    echo "Running all build"
    ${SCRIPT_DIR}/make_all.sh $1 $2
}

function fluya_build(){
    echo "Running fluya build"
    ${SCRIPT_DIR}/make_fluya.sh $1 $2
}

function local_build(){
    echo "Running local build"
    ${SCRIPT_DIR}/make_local_install.sh $1 $2
}

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments provided, try using ./utopia_compile help."
    exit 1
fi

arg2="1"
arg3="/usr/local/"

if [[ "$OSTYPE" == "win32" ]]; then
    arg3="c:/Program Files/utopia"
fi




while getopts ":b:j:p:h" opt; do
  case $opt in
    b) arg1="$OPTARG"
    ;;
    j) arg2="$OPTARG"
    ;;
    h) display_help
        exit 0
    ;;
    p) arg3="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done


case $arg1 in
    basic)
        basic_build $arg2 $arg3
        ;;
    all)
        all_build $arg2 $arg3
        ;;
    fluya)
        fluya_build $arg2 $arg3
        ;;
    local)
        local_build $arg2 $arg3
        ;;
    *)
        echo "Unknown build type: $build_type"
        ;;
esac

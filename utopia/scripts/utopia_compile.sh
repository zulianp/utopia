#!/bin/bash

# Function to display help text
function display_help() {
    echo "This is the help text."
}

# Function for basic build
function basic_build() {
    echo "Running basic build..."
    echo "${BASH_SOURCE}"
    .$PWD/make_basic.sh
}

# Function for all build
function all_build() {
    echo "Running all build..."
    .$PWD/make_all.sh
}

# Check if an argument was provided
if [ $# -eq 0 ]; then
    echo "No arguments provided. Use './example.sh help' for help."
    exit 1
fi

# Parse arguments
for arg in "$@"
do
    case $arg in
        help)
            display_help
            ;;
        build=*)
            build_type="${arg#*=}"
            case $build_type in
                basic)
                    basic_build
                    ;;
                all)
                    all_build
                    ;;
                *)
                    echo "Unknown build type: $build_type"
                    ;;
            esac
            ;;
        *)
            echo "Unknown argument: $arg"
            ;;
    esac
done

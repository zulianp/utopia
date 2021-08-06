#!/bin/bash

if [ -n "$INSTALL_DIR" ]; then
    echo "INSTALL_DIR="$INSTALL_DIR
else
    echo "define INSTALL_DIR where you want trilinos to be installed"
    return
fi


if [ -n "$UTOPIA_SRC_DIR" ]; then
    echo "UTOPIA_SRC_DIR="$UTOPIA_SRC_DIR
else
    echo "define UTOPIA_SRC_DIR to the root of utopia git"
    return
fi


git clone https://github.com/trilinos/Trilinos.git; \
    cd Trilinos && \
    mkdir build; \
    cd build && \
    source $UTOPIA_SRC_DIR/utopia/scripts/installer/configure_trilinos.sh \
    make -j && make install
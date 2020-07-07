#!/bin/bash

# Copy me into your bin folder

# Uncomment me to remove previous coverage runs
rm -f ./*.info; find . -name "*.gcda" -print0 | xargs -0 rm

lcov --capture --directory . --output-file cov.info

# Make sure we have right configuration
cmake .. -DCMAKE_BUILD_TYPE=Debug -DENABLE_DEPRECATED_API=OFF
make -j4 complete && ./utopia_test -verbose &&  ./utopia_test -verbose -test unconstrained_opt -test newton_ls -test constrained_opt && ./utopia_bench -verbose
ret=$?

echo 'execution returned' $ret

if [[ $ret -eq '0' ]]; then
    echo "Tests successful"
else
    echo "Did not compile or pass the test";
    return 1;
fi

lcov --capture --directory . --output-file cov.info
lcov --remove cov.info '*/utopia/apps/*' '/usr/include/*' '/usr/lib/*' '/Applications/*' '*/utopia/external/*' '*/rapidxml/*' '*/include/*' -o filtered_cov.info
genhtml  --prefix ../ --ignore-errors source filtered_cov.info  --output-directory out

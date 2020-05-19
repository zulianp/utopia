#!/bin/bash

# Copy me into your bin folder

# Uncomment me to remove previous coverage runs
rm ./*.info; find . -name "*.gcda" -print0 | xargs -0 rm

lcov --capture --directory . --output-file cov.info

make -j4 complete && ./utopia_test -verbose && ./utopia_bench
ret=$?

echo 'execution returned' $ret

if [[ $ret -eq '0' ]]; then
    echo "Tests successful"
else
    echo "Did not compile or pass the test";
    return;
fi

lcov --capture --directory . --output-file cov.info
lcov --remove cov.info '*/utopia/apps/*' '/usr/include/*' '/usr/lib/*' '/Applications/*' '*/utopia/external/*' '*/rapidxml/*' '*/include/*' -o filtered_cov.info
genhtml  --prefix ../ --ignore-errors source filtered_cov.info  --output-directory out

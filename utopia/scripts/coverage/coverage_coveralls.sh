#!/bin/bash
coveralls --exclude /usr --gcov-options '\-lp' --root $PWD/..;
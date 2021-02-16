#!/bin/bash

docker build -t utopia .
docker run -it -v /Users/cesana/desktop/tesi:/utopia/utopia/build/scripting/tesi utopia sh


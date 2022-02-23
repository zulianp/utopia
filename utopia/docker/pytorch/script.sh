#!/bin/bash

docker build -t utopia .
docker run -it -v /Users/cesana/desktop/utopia/utopia:/work_space utopia sh


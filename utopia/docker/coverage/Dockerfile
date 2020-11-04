FROM zulianp/utopia-env-complete-edge:latest

# Update to last version
RUN cd /utopia_cloned/utopia && git checkout development && git pull
RUN rm -rf /utopia_cloned/utopia/bin_cov

# Compile
RUN cd /utopia_cloned/utopia && mkdir bin_cov
RUN cd /utopia_cloned/utopia/bin_cov && cmake -DUTOPIA_ENABLE_CODE_COVERAGE=ON -DCMAKE_BUILD_TYPE=Debug ..
RUN cd /utopia_cloned/utopia/bin_cov && make complete

# Create coverage data
RUN cp /utopia_cloned/utopia/scripts/utils/coverage.sh /utopia_cloned/utopia/bin_cov/coverage.sh
RUN cd /utopia_cloned/utopia/bin_cov && source coverage.sh

# Running in parallel inside Docker
RUN echo "For parallel runs use: mpirun --allow-run-as-root -n 2 ./utopia_test"


RUN cd /utopia_cloned/utopia/bin_cov && CODECOV_TOKEN=c754a0a2-6280-4ce4-9756-b65d1f910144 bash <(curl -s https://codecov.io/bash)
CMD ["/bin/sh"]

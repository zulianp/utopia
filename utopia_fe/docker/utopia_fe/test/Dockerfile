ARG DOCKER_TAG

FROM zulianp/utopia:$DOCKER_TAG
RUN echo 'tag='${DOCKER_TAG}

# On laptop
# FROM zulianp/utopia:latest
# docker-compose -f docker-compose.test.yml run sut

COPY . .
RUN ls

CMD ["/bin/sh"]
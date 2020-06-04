
# ARG IMAGE_NAME
# FROM ${IMAGE_NAME}
# ARG IMAGE_NAME

# This only works for Pull requests
FROM zulianp/utopia:this

RUN echo 'IMAGE_NAME='${IMAGE_NAME}

# On laptop
# FROM zulianp/utopia:latest
# docker-compose -f docker-compose.test.yml run sut

COPY . .
RUN ls

CMD ["/bin/sh"]
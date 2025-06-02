FROM n8nio/n8n:latest

USER root

RUN apk update && \
    apk add docker-cli

RUN addgroup -g 972 docker && adduser node docker
USER node

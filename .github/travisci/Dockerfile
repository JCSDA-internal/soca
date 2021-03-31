FROM  jcsda/docker-gnu-openmpi-dev:latest

RUN  mkdir -p /jcsda \
    && chown -R jedi:jedi /jcsda \
    && chmod 6755 /jcsda

USER jedi

COPY default-mca-params.conf /jcsda/mca-params.conf

RUN mkdir -p ~/.openmpi \
    && cp /jcsda/mca-params.conf ~/.openmpi/

CMD ["/bin/bash"]

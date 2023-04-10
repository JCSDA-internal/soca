FROM  jcsda/docker-gnu-openmpi-dev:latest

RUN  mkdir -p /jcsda \
    && chown -R nonroot:nonroot /jcsda \
    && chmod 6755 /jcsda

USER nonroot

COPY default-mca-params.conf /jcsda/mca-params.conf

RUN mkdir -p ~/.openmpi \
    && cp /jcsda/mca-params.conf ~/.openmpi/

CMD ["/bin/bash"]

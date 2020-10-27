FROM  jcsda/docker-gnu-openmpi-dev:latest

RUN groupadd jcsda -g 9999 \
    && adduser jcsdauser \
    && mkdir -p /jcsda \
    && chown -R jcsdauser:jcsda /jcsda \
    && chmod 6755 /jcsda

USER jcsdauser

COPY default-mca-params.conf /jcsda/mca-params.conf

RUN mkdir -p ~/.openmpi \
    && cp /jcsda/mca-params.conf ~/.openmpi/

CMD ["/bin/bash"]

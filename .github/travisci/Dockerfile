FROM  jcsda/docker-gnu-openmpi-dev:latest

RUN cd /usr/local/src \
    && curl -L -O https://github.com/ccache/ccache/releases/download/v3.7.4/ccache-3.7.4.tar.gz \
    && tar -xaf ccache-3.7.4.tar.gz \
    && cd ccache-3.7.4 \
    && ./configure \
    && make \
    && make install

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

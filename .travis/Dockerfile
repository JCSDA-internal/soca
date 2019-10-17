FROM  jcsda/docker:latest

RUN cd /usr/local/src \
    && curl -L -O https://github.com/ccache/ccache/releases/download/v3.7.4/ccache-3.7.4.tar.gz \
    && tar -xaf ccache-3.7.4.tar.gz \
    && cd ccache-3.7.4 \
    && ./configure \
    && make \
    && make install
    
RUN mkdir -p /var/run/sshd \
    && ssh-keygen -A \
    && sed -i 's/#PermitRootLogin yes/PermitRootLogin yes/g' /etc/ssh/sshd_config \
    && sed -i 's/#RSAAuthentication yes/RSAAuthentication yes/g' /etc/ssh/sshd_config \
    && sed -i 's/#PubkeyAuthentication yes/PubkeyAuthentication yes/g' /etc/ssh/sshd_config

RUN groupadd jcsda -g 9999 \
    && adduser jcsdauser \
    && mkdir -p /jcsda \
    && chown -R jcsdauser:jcsda /jcsda \
    && chmod 6755 /jcsda 

RUN mkdir /jcsda/.ssh \
    && echo "StrictHostKeyChecking no" > /jcsda/.ssh/config 
RUN chown -R jcsdauser:jcsda /jcsda/

USER jcsdauser
WORKDIR /jcsda

RUN ssh-keygen -f /jcsda/.ssh/id_rsa -t rsa -N '' \
    && chmod 600 /jcsda/.ssh/config \
    && chmod 700 /jcsda/.ssh \
    && cp /jcsda/.ssh/id_rsa.pub /jcsda/.ssh/authorized_keys

VOLUME /jcsda

CMD ["/bin/bash"]

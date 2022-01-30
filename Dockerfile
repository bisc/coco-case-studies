FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    bsdmainutils \
    less \
    locales \
    locales-all \
    unzip \
    wget \
    xz-utils \      
    bsdmainutils \
    less \
    vim \
    nano \
    && rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

RUN wget https://account.wolfram.com/download/public/wolfram-engine/desktop/LINUX && bash LINUX -- -auto -verbose && rm LINUX

RUN echo "function pretty_csv { \n\tcolumn -t -s, -n \"\$@\" | less -F -S -X -K \n}" >> /root/.bashrc

COPY coco-case-studies /coco-case-studies
WORKDIR /coco-case-studies

COPY ./docker-entrypoint.sh /
ENTRYPOINT ["/docker-entrypoint.sh"]
CMD ["/bin/bash"]

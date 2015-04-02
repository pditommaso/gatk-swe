FROM pditommaso/dkrbase
MAINTAINER Paolo Di Tommaso

RUN apt-get update --fix-missing && \
    apt-get install -y bzip2 libncurses5-dev

RUN wget -q 'http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbio-bwa%2Ffiles%2F&ts=1427965175&use_mirror=garr' -O - | tar xj && \
 cd bwa-0.7.12 && \
 make && \
 mv bwa /usr/local/bin/

RUN wget -q -O- 'http://downloads.sourceforge.net/project/samtools/samtools/1.1/samtools-1.1.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2F&ts=1420911249&use_mirror=heanet' | tar xj  && \
 cd samtools-1.1 && \
 make && make install && \
 cd htslib-1.1/ && \
 make && make install

RUN wget -q http://genome.crg.es/~epalumbo/gene2farm/GenomeAnalysisTK-2.7-4.tar.bz2 -O- | tar xj

RUN apt-get install -y openjdk-7-jre-headless
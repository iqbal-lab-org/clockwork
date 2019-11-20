
FROM ubuntu:18.04

## need the following environment variable for automated installation of tzdata
ENV DEBIAN_FRONTEND=noninteractive 


LABEL authors="Martin Hunt (ebi.ac.uk), Evan Biederstedt (evan.biederstedt@gmail.com)" \
    version.image="0.1.0" \
    version.clockwork="v0.7.8" \
    version.bcftools="1.3.1" \
    version.bwa="v0.7.15" \
    version.perl="5.24.1" \
    version.enaBrowserTools="v1.4.1" \
    version.fastqc="v0.11.5." \
    version.picard="2.9.4" \
    version.seqtk="v1.2" \
    version.samtools="1.3.1" \
    version.trimmomatic="0.36" 


## from https://github.com/iqbal-lab-org/clockwork/blob/master/scripts/install_dependencies.sh

# NOTE: needed to specify java version 8 (openjdk-8-jre) because
# default is version 10 (from default-jre), which won't work with nextflow.


## set Perl environment variables
ENV PATH="/bioinf-tools/:/bioinf-tools/cortex/scripts/analyse_variants/seq-align/bin/:/bioinf-tools/vcftools-0.1.15/install/bin:/bioinf-tools/enaBrowserTools-1.4.1/python3:/clockwork/scripts/:$PATH"
ENV PERL5LIB="/bioinf-tools/vcftools-0.1.15/install/share/perl/5.24.1/:/bioinf-tools/cortex/scripts/analyse_variants/bioinf-perl/lib:/bioinf-tools/cortex/scripts/calling/:$PERL5LIB"
ENV LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"



RUN apt-get update
RUN apt-get install -y software-properties-common
RUN apt-add-repository universe
RUN apt-get update

RUN apt-get install -y \
    build-essential \
    cmake \
    curl \
    gawk \
    git \
    gnuplot \
    graphviz \
    openjdk-8-jre \
    libarchive-dev \
    liblzma-dev \
    libbz2-dev \
    libhts-dev \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    pkg-config \
    python-dev \
    python-pip \
    python3-dev \
    python3-pip \
    python3-setuptools \
    r-base-core \
    rsync \
    unzip \
    wget --fix-missing


#_________________________ bcftools _______________________#

RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 \
    && tar xf bcftools-1.3.1.tar.bz2 \
    && cd bcftools-1.3.1/ \
    && make \ 
    && cd .. \
    && cp bcftools-1.3.1/bcftools /bin/bcftools


#__________________________ BWA____________________________#
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 \
    && tar xf bwa-0.7.15.tar.bz2 \
    && cd bwa-0.7.15/ \
    && make \
    && cd .. \
    && cp bwa-0.7.15/bwa /bin/bwa


#_____________________ enaBrowserTools ____________________#
RUN wget https://github.com/enasequence/enaBrowserTools/archive/v1.4.1.tar.gz \
    && tar xf v1.4.1.tar.gz


#_________________________ FASTQC ________________________#
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip \
    && unzip fastqc_v0.11.5.zip \
    && chmod 755 FastQC/fastqc \
    && cp FastQC/fastqc /bin/fastqc

#_________________________ qctools ________________________#
RUN wget https://github.com/alastair-droop/fqtools/archive/986e451.tar.gz \
    && tar xf 986e451.tar.gz \
    && rm 986e451.tar.gz \
    && cd fqtools-986e451/ \
    && make \
    && cd .. \
    && cp fqtools-986e451/bin/fqtools /bin/fqtools

#________________________ mccortex _______________________#
RUN git clone --recursive https://github.com/mcveanlab/mccortex \
    && cd mccortex \
    && git checkout 97aba198d632ee98ac1aa496db33d1a7a8cb7e51 \
    && make all \
    && cd .. \
    && cp mccortex/bin/mccortex31 /bin/mccortex31 


#________________________ Mykrobe ________________________#
RUN git clone https://github.com/Mykrobe-tools/mykrobe.git mykrobe \
    && cd mykrobe \ 
    && git checkout b47f6241cfb2d223bb53aafd90a3f91d97248052  \ 
    && sed -i 's/help="build variant probes", aliases=.*$/help="build variant probes")/' src/mykrobe/cli.py \
     # fix for python2: aliases is not an option so remove it
    && wget -O mykrobe-data.tar.gz https://goo.gl/DXb9hN \ 
    && tar -zxvf mykrobe-data.tar.gz && rm -fr src/mykrobe/data \ 
    && mv mykrobe-data src/mykrobe/data \
    && pip3 install .


#________________________ nextflow _______________________#
RUN wget -qO- get.nextflow.io | bash \
    && chmod 755 nextflow

#________________________ picard _________________________#
RUN wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar

#________________________ seqtk __________________________#
RUN wget https://github.com/lh3/seqtk/archive/v1.2.tar.gz \
    && tar xf v1.2.tar.gz \
    && rm v1.2.tar.gz \
    && cd seqtk-1.2/ \
    && make \
    && cd .. \
    && cp seqtk-1.2/seqtk /bin/seqtk

#_________________________ samtools ______________________#
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 \
    && tar xf samtools-1.3.1.tar.bz2 \
    && cd samtools-1.3.1/ \
    && make \
    && cd .. \
    && cp samtools-1.3.1/samtools /bin/samtools \
    && cp -rp samtools-1.3.1/misc/plot-bamstats /bin/plot-bamstats

#________________________ stampy _________________________#
RUN wget http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz \
    && tar xf stampy-latest.tgz \
    && rm stampy-latest.tgz \
    && cd stampy-* \
    && make \
    && cd .. \
    && cp stampy-*/stampy.py /bin/stampy.py

#________________________ Trimmomatic ____________________#
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip \
    && unzip Trimmomatic-0.36.zip

#________________________ vcftools _______________________#
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz \
    && tar xf vcftools-0.1.15.tar.gz \
    && cd vcftools-0.1.15 \
    && ./configure --prefix $PWD/install \
    && make \
    && make install

# cortex needs the perl/ directory. It expects it to be in the vcftools root,
# but somehwere between v0.1.9 and v0.1.15 it moved into src/.
RUN ln -s src/perl/ .

#________________________ cortex _________________________#
RUN wget --no-check-certificate -O cortex.tar.gz https://github.com/iqbal-lab/cortex/archive/master.tar.gz \
    && tar xf cortex.tar.gz \
    && mv cortex-master cortex \
    && cd cortex/ \
    && bash install.sh \
    && make NUM_COLS=1 cortex_var \
    && make NUM_COLS=2 cortex_var

# ___________________ python packages ___________________#
# note: requests needs to be here instead of as part of
# python setup.py install, because setup.py install
# throws an error.  This way works.
RUN pip3 install python-dateutil requests pysam pyfastaq pymysql numpy openpyxl pyflakes scipy XlsxWriter



#________________________ gramtools _________________________#
RUN pip3 install git+https://github.com/iqbal-lab-org/gramtools@4644135e162ff1ff333f970f1dafd919aa0fc279

#________________________ mummer ____________________________#
RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz  \
    && tar xf mummer-4.0.0beta2.tar.gz \
    && cd mummer-4.0.0beta2 \
    && ./configure \
    && make \
    && make install


#________________________ minos _____________________________#
RUN pip3 install bio-minos==0.9.1






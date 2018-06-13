#!/usr/bin/env bash
set -vexu

install_root=$1


# We need packages from universe, so make sure it's added
# (it isn't for singularity)
apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update
apt-get install -y \
  build-essential \
  cmake \
  curl \
  default-jre \
  gawk \
  git \
  gnuplot \
  graphviz \
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
  wget


mkdir $install_root
cd $install_root

#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2
tar xf bcftools-1.3.1.tar.bz2
cd bcftools-1.3.1/
make
cd ..
cp -s bcftools-1.3.1/bcftools .


#__________________________ BWA____________________________#
cd $install_root
wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2
tar xf bwa-0.7.15.tar.bz2
cd bwa-0.7.15/
make
cd ..
cp -s bwa-0.7.15/bwa .


#_____________________ enaBrowserTools ____________________#
cd $install_root
wget https://github.com/enasequence/enaBrowserTools/archive/v1.4.1.tar.gz
tar xf v1.4.1.tar.gz


#_________________________ FASTQC ________________________#
cd $install_root
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
cp -s FastQC/fastqc .

#_________________________ qctools ________________________#
cd $install_root
wget https://github.com/alastair-droop/fqtools/archive/986e451.tar.gz
tar xf 986e451.tar.gz
rm 986e451.tar.gz
cd fqtools-986e451/
make
cd ..
cp -s fqtools-986e451/bin/fqtools .

#________________________ mccortex _______________________#
cd $install_root
git clone --recursive https://github.com/mcveanlab/mccortex
cd mccortex
git checkout 97aba198d632ee98ac1aa496db33d1a7a8cb7e51
make all
cd ..
cp -s mccortex/bin/mccortex31 .


#________________________ Mykrobe ________________________#
cd $install_root
git clone https://github.com/Mykrobe-tools/mykrobe-atlas-cli.git mykrobe
cd mykrobe
#git checkout fd3759a1f30aec61346b7fd7888df8e8e9d3025e

# This is new branch with options to use own panels/probes/json:
git checkout e96cdb0815575dc9fb1c91b20235d048dd5ce4f0
# fix for python2: aliases is not an option so remove it
sed -i 's/help="build variant probes", aliases=.*$/help="build variant probes")/' src/mykrobe/cli.py
wget -O mykrobe-data.tar.gz https://goo.gl/DXb9hN && tar -zxvf mykrobe-data.tar.gz && rm -fr src/mykrobe/data && mv mykrobe-data src/mykrobe/data
pip install .


#________________________ nextflow _______________________#
cd $install_root
wget -qO- get.nextflow.io | bash
chmod 755 nextflow

#________________________ picard _________________________#
cd $install_root
wget https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar

#________________________ seqtk __________________________#
cd $install_root
wget https://github.com/lh3/seqtk/archive/v1.2.tar.gz
tar xf v1.2.tar.gz
rm v1.2.tar.gz
cd seqtk-1.2/
make
cd ..
cp -s seqtk-1.2/seqtk .

#_________________________ samtools ______________________#
cd $install_root
wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2
tar xf samtools-1.3.1.tar.bz2
cd samtools-1.3.1/
make
cd ..
cp -s samtools-1.3.1/samtools .
cp -rp samtools-1.3.1/misc/plot-bamstats .

#________________________ stampy _________________________#
cd $install_root
wget http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz
tar xf stampy-latest.tgz
rm stampy-latest.tgz
cd stampy-*
make
cd ..
cp -s stampy-*/stampy.py .

#________________________ Trimmomatic ____________________#
cd $install_root
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip

#________________________ vcftools _______________________#
cd $install_root
wget https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz
tar xf vcftools-0.1.15.tar.gz
cd vcftools-0.1.15
./configure --prefix $PWD/install
make
make install

# cortex needs the perl/ directory. It expects it to be in the vcftools root,
# but somehwere between v0.1.9 and v0.1.15 it moved into src/.
ln -s src/perl/ .

#________________________ cortex _________________________#
cd $install_root
wget --no-check-certificate -O cortex.tar.gz https://github.com/iqbal-lab/cortex/archive/master.tar.gz
tar xf cortex.tar.gz
mv cortex-master cortex
cd cortex/
bash install.sh
make NUM_COLS=1 cortex_var
make NUM_COLS=2 cortex_var

# ___________________ python packages ___________________#
# note: requests needs to be here instead of as part of
# python setup.py install, because setup.py install
# throws an error.  This way works.
pip3 install python-dateutil requests pysam pyfastaq pymysql numpy openpyxl pyflakes scipy XlsxWriter



#________________________ gramtools _________________________#
pip3 install git+https://github.com/iqbal-lab-org/gramtools@03ff4ff2309f25f5470816b3097d60ea7f0a9754

#________________________ mummer ____________________________#
cd $install_root
wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar xf mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure
make
make install


#________________________ minos _____________________________#
pip3 install git+https://github.com/iqbal-lab-org/minos@9662d64ef83b67b62372b569d994eea28e7cbb63


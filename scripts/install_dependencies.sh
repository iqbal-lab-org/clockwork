#!/usr/bin/env bash
set -vexu

install_root=$1


# We need packages from universe, so make sure it's added
# (it isn't for singularity)
apt-get update
apt-get install -y software-properties-common
apt-add-repository universe
apt-add-repository multiverse
apt-get update
apt-get upgrade
echo "__________________ R version __________________________"
R --version
DEBIAN_FRONTEND=noninteractive apt-get install -y \
  automake \
  build-essential \
  cmake \
  curl \
  gawk \
  git \
  gnuplot \
  graphviz \
  openjdk-8-jre \
  libarchive-dev \
  libcurl4-gnutls-dev \
  liblzma-dev \
  libbz2-dev \
  libhts-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libvcflib-tools \
  libssl-dev \
  zlib1g-dev \
  pkg-config \
  python3-dev \
  python3-pip \
  python3-setuptools \
  rsync \
  unzip \
  tabix \
  wget

# Note: needed to specify java version 8 (openjdk-8-jre) because
# default is version 10 (from default-jre), which won't work with nextflow.

mkdir $install_root
cd $install_root

#_________________________ bcftools _______________________#
cd $install_root
wget -q https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2
tar xf bcftools-1.15.1.tar.bz2
cd bcftools-1.15.1/
make
cd ..
cp -s bcftools-1.15.1/bcftools .

#_____________________ enaBrowserTools ____________________#
cd $install_root
git clone https://github.com/enasequence/enaBrowserTools.git
cd enaBrowserTools
git checkout 7075a896f822e3ea3d3fac8bc10bcfeeb2506685


#_________________________ FASTQC ________________________#
cd $install_root
wget -q https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 FastQC/fastqc
cp -s FastQC/fastqc .

#_________________________ qctools ________________________#
cd $install_root
wget -q https://github.com/alastair-droop/fqtools/archive/986e451.tar.gz
tar xf 986e451.tar.gz
rm 986e451.tar.gz
cd fqtools-986e451/
make
cd ..
cp -s fqtools-986e451/bin/fqtools .

#________________________ mccortex _______________________#
cd $install_root
git clone --recursive https://github.com/Mykrobe-tools/mccortex
cd mccortex
git checkout 5a9d410468f6b2980434e415ec341be320d37d82
make all
cd ..
cp -s mccortex/bin/mccortex31 .

#------------------------------ minimap2 ---------------------------------------
cd $install_root
MINIMAP2_V=2.24
wget https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_V}/minimap2-${MINIMAP2_V}.tar.bz2
tar xf minimap2-${MINIMAP2_V}.tar.bz2
rm minimap2-${MINIMAP2_V}.tar.bz2
cd minimap2-${MINIMAP2_V}
arch_is_arm=$(dpkg --print-architecture | grep '^arm' | wc -l)
if [[ $arch_is_arm -gt 0 ]]
then
    make arm_neon=1 aarch64=1
else
    make
fi
cd ..
cp -s minimap2-${MINIMAP2_V}/minimap2 .

#________________________ Mykrobe ________________________#
cd $install_root
git clone https://github.com/Mykrobe-tools/mykrobe.git mykrobe
cd mykrobe
git checkout 17540f25d6b84b5cb3dfa59973d5838f1e1cb51c
pip3 install requests
rm -rf mccortex
git clone --recursive -b geno_kmer_count https://github.com/Mykrobe-tools/mccortex mccortex
cd mccortex
make
cd ..
pip3 install .
myk_dir=$(pip3 show mykrobe | awk '/^Location/ {print $NF}')
echo $myk_dir
cp mccortex/bin/mccortex31 $myk_dir/mykrobe/cortex/mccortex31
mykrobe panels update_metadata --debug
mykrobe panels update_species --debug all


#________________________ nextflow _______________________#
cd $install_root
wget -qO- get.nextflow.io | bash
chmod 755 nextflow

#________________________ picard _________________________#
cd $install_root
wget -q https://github.com/broadinstitute/picard/releases/download/2.9.4/picard.jar

#________________________ seqtk __________________________#
cd $install_root
wget -q https://github.com/lh3/seqtk/archive/v1.2.tar.gz
tar xf v1.2.tar.gz
rm v1.2.tar.gz
cd seqtk-1.2/
make
cd ..
cp -s seqtk-1.2/seqtk .

#_________________________ samtools ______________________#
cd $install_root
wget -q https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2
tar xf samtools-1.15.1.tar.bz2
cd samtools-1.15.1/
make
cd ..
cp -s samtools-1.15.1/samtools .
cp -rp samtools-1.15.1/misc/plot-bamstats .


#________________________ Trimmomatic ____________________#
cd $install_root
wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip

#________________________ vcftools _______________________#
cd $install_root
wget -q https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz
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
git clone --recursive https://github.com/iqbal-lab/cortex.git
cd cortex
git checkout c8147152cd4015c45057900e8fb600376d1d7fb3
bash install.sh
make NUM_COLS=1 cortex_var
make NUM_COLS=2 cortex_var

# ___________________ python packages ___________________#
# note: requests needs to be here instead of as part of
# python setup.py install, because setup.py install
# throws an error.  This way works.
pip3 install cython
pip3 install python-dateutil requests pysam pyfastaq pymysql numpy openpyxl pyflakes scipy XlsxWriter



#________________________ gramtools _________________________#
cd $install_root
# Why six>=1.14.0?
# See https://github.com/pypa/virtualenv/issues/1551
pip3 install tox "six>=1.14.0"
git clone https://github.com/iqbal-lab-org/gramtools
cd gramtools
git checkout 8af53f6c8c0d72ef95223e89ab82119b717044f2
# Note: a simple "pip3 install ." works for singularity but
# not for docker - the `gram` exectuable does not get
# put where gramtools expects to find it. The method
# below, which explicitly builds the binary, then installs
# does work ok for both docker and singularity.
mkdir cmake-build
cd cmake-build
cmake .. -DCMAKE_BUILD_TYPE=REL_WITH_ASSERTS
make gram
cd ..
pip3 install -e .

#________________________ mummer ____________________________#
cd $install_root
wget -q https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz
tar xf mummer-4.0.0beta2.tar.gz
cd mummer-4.0.0beta2
./configure
make
make install


#________________________ vt __________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .

#______________________ ivcmerge ______________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/ivcfmerge.git
cd ivcfmerge
git checkout 5819787614a263a9f35fd0c247442f092ab174ff
pip3 install .

#________________________ minos _____________________________#
pip3 install 'cluster_vcf_records==0.13.2'
pip3 install git+https://github.com/iqbal-lab-org/minos@v0.12.0

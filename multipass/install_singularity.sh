#!/usr/bin/env bash

apt install -y \
    autoconf \
    build-essential \
    libssl-dev \
    uuid-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    pkg-config

arch_is_arm=$(dpkg --print-architecture | grep '^arm' | wc -l)

# ------------------------------------- go -------------------------------------
GO_VERSION=1.18.1
if [[ $arch_is_arm -gt 0 ]]
then
    GO_TAR=go${GO_VERSION}.linux-arm64.tar.gz
else
    GO_TAR=go${GO_VERSION}.linux-amd64.tar.gz
fi
wget https://go.dev/dl/${GO_TAR}
rm -rf /usr/local/go
tar -C /usr/local -xzf ${GO_TAR}
rm ${GO_TAR}
export PATH=$PATH:/usr/local/go/bin


# --------------------------------- singularity --------------------------------
export SING_VERSION=3.9.9
export GOPATH=/tmp/go
mkdir -p $GOPATH/src/github.com/sylabs
cd $GOPATH/src/github.com/sylabs
wget https://github.com/sylabs/singularity/releases/download/v${SING_VERSION}/singularity-ce-${SING_VERSION}.tar.gz
tar -xzf singularity-ce-${SING_VERSION}.tar.gz
cd singularity-ce-${SING_VERSION}
./mconfig
make -C ./builddir
make -C ./builddir install
cd $HOME
rm -rf $GOPATH

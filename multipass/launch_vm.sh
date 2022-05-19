#!/usr/bin/env bash

if [ $# -ne 1 ]
then
    echo "usage: $0 name_of_vm"
    echo
    exit
fi

name=$1

multipass launch -m 10G -d 20G -n ${name} 20.04
multipass mount $HOME ${name}:/home/ubuntu/Home

multipass exec ${name} sudo apt-get update
multipass exec ${name} -- sudo apt-get upgrade -y

multipass transfer install_mysql.sh ${name}:.
multipass exec ${name} sudo bash install_mysql.sh

multipass transfer install_singularity.sh ${name}:.
multipass exec ${name} sudo bash install_singularity.sh

multipass transfer ../scripts/install_dependencies.sh ${name}:.
multipass exec ${name} sudo bash install_dependencies.sh /bioinf-tools

multipass transfer bashrc ${name}:.bashrc_user
multipass exec ${name} -- bash -c "echo source \$HOME/.bashrc_user >> ~/.bashrc"

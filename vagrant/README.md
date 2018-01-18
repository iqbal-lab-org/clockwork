# Clockwork Vagrant Virtual Machine

This is for developers to create a VM to develop the code and run the
tests. Assumes you have vagrant installed.

Note: you will probably want to change the lines in `VagrantFile` that mount
directories from the host inside the VM.

Set up, start, and log into the VM by running these
commands from this (`vagrant/`) directory:

    vagrant up
    vagrant ssh

Then inside the VM, clone this repository, change to the `python/`
directory in the repository and run the tests with:

    python3 setup.py test


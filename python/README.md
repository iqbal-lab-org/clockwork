# Clockwork Python code

This is all intended to be run inside the clockwork singularity container.


If you want to develop the code and/or run the tests, please look
in the `vagrant/` directory for how to make a VM to do this.

## I really want to install everything by hand anyway

You will need these dependencies:

* BWA
* bcftools
* cortex. Set the environment variable `CLOCKWORK_CORTEX_DIR` to
  the absolute path to the root cortex directory.
* enaBrowserTools (so the scripts enaDataGet and enaGroupGet work)
* fastqc
* fqtools
* picard. Set the environment variable `CLOCKWORK_PICARD_JAR` to
  the absolute path to the picard .jar file.
* samtools
* stampy. Set the environment variable `CLOCKWORK_STAMPY_SCRIPT` to
  the absoulate path to the `stampy.py` script.
* trimmomatic. Set the environment variable `CLOCKWORK_TRIMMO_DIR` to
  the absolute path to the root Trimmomatic directory.
* vcftools. Set the environment variable `CLOCKWORK_VCFTOOLS_DIR` to
  the absolute path to the root vcftools directory.


The code will check first for each of those environment variables,
Whenever one is not defined, the code then looks for
`/bioinf-tools/` for executables/directories/jar files.
The vagrant VM and the container
install everything in `/bioinf-tools`. The script
`scripts/install_dependencies.sh` may help - it is used to install
dependencies in the Vagrant VM and in the Singularity container.



Run the tests:

    python3 setup.py test

If the tests pass, install:

    python3 setup.py install



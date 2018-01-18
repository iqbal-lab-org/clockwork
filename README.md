# Clockwork - pipelines for processing bacteria Illumina data and variant calling

Note: these pipelines were made for the [CRyPTIC](http://www.crypticproject.org/) project, but
in principle can be used on any bacteria.

## Installation

### Dependencies

* [Singularity](http://singularity.lbl.gov/)
* [MySQL](https://www.mysql.com/) (not needed if you want to run pipelines on each sample individually)
* [nextflow](https://www.nextflow.io/)


### Build the Singularity container

Run this to make the container `clockwork_container.img` (you must be inside the
directory `singularity/` for the build to work):

    cd singularity/
    sudo singularity create -s 8000 clockwork_container.img
    sudo singularity bootstrap clockwork_container.img clockwork_container.def


### Make pipeline directories (only if using a database)

Make the following directories (you can call them what you like,
but we will use these names throughout this document):

    mkdir Pipeline_root # this is where reads and analysis files are kept long-term
    mkdir Pipeline_refs # where reference genomes are stored
    mkdir Pipeline_spreadsheet_archive # import spreadsheets archived in here

### Set up database

You need to make a config file with the login and database name. Here is an example:

```
[db_login]
user = user_name
password = password
host = host_name
db = name_of_database
```

There can optionally be a line for the port: `port = N`. If this is not
present then the default port is used.

Assuming your config file is called `db.ini`, run this to set up the tracking
database:

    singularity exec clockwork_container.img clockwork make_empty_db db.ini

### Set up reference genomes

The pipelines require the following reference sequences:

1. Genomes to use for removing contamination from reads.
2. At least one reference genome to use for QC and variant calling.

The container has references for TB, or you can provide your own (to
be documented).

To set up the contamination reference if you are using a database:

    singularity exec clockwork_container.img clockwork reference_prepare \
      --db_config_file db.ini --pipeline_references_root Pipeline_refs \
      --contam_tsv /clockwork/reference_genomes/remove_contam.tsv \
      --name remove_contam \
      /clockwork/reference_genomes/remove_contam.fa.gz

Note: the files in `/clockwork/` are in the container. You should not
expect them to be present on your host.

Or, if you are not using a database:

    singularity exec clockwork_container.img clockwork reference_prepare \
      --contam_tsv /clockwork/reference_genomes/remove_contam.tsv \
      --outdir Reference.remove_contam \
      /clockwork/reference_genomes/remove_contam.fa.gz

To setup the QC/variant reference using a database:

    singularity exec clockwork_container.img clockwork reference_prepare \
      --db_config_file db.ini --pipeline_references_root Pipeline_refs \
      --name NC_000962.3 \
      /clockwork/reference_genomes/NC_000962.3.fa


Or, if you are not using a database:

    singularity exec clockwork_container.img clockwork reference_prepare \
      --outdir Reference.NC_000962.3 \
      /clockwork/reference_genomes/NC_000962.3.fa


## Import pipeline

This is only relevant if you are using a database. It imports FASTQ files,
taking metadata from a spreadsheet (in xlsx or TSV format). Run:

    nextflow run nextflow/import.nf \
      -with-singularity clockwork_container.img \
      --dropbox_dir dir_with_spreadheet_and_fastq_files \
      --pipeline_root Pipeline_root \
      --xlsx_archive_dir Pipeline_spreadsheet_archive \
      --db_config_file db.ini \

It will find all spreadsheets called `*.xlsx` or `*.tsv`, and import their
data.

## Remove contamination pipeline

### Using a database

    nextflow run nextflow/remove_contam.nf \
      -with-singularity clockwork_container.img \
      --ref_id 1 \
      --references_root Pipeline_refs \
      --pipeline_root Pipeline_root \
      --db_config_file db.ini

### Run on one read pair

    nextflow run nextflow/remove_contam.nf \
      -with-singularity clockwork_container.img \
      --ref_fasta Reference.remove_contam/ref.fa \
      --ref_metadata_tsv Reference.remove_contam/remove_contam_metadata.tsv \
      --reads_in1 reads_1.fastq --reads_in2 reads_2.fastq \
      --outprefix OUT


## QC pipeline

### Using a database

    nextflow run nextflow/qc.nf \
      -with-singularity clockwork_container.img \
      --ref_id 2 \
      --references_root  Pipeline_refs \
      --pipeline_root Pipeline_root \
      --db_config_file db.ini

### Run on one read pair

    nextflow run nextflow/qc.nf \
      -with-singularity clockwork_container.img \
      --ref_fasta Reference.NC_000962.3/ref.fa \
      --reads_in1 reads_1.fastq --reads_in2 reads_2.fastq \
      --output_dir qc_out


## Variant calling pipeline

### Using a database

    nextflow run nextflow/variant_call.nf \
      -with-singularity clockwork_container.img \
      --ref_id 2 \
      --references_root  Pipeline_refs \
      --pipeline_root Pipeline_root \
      --db_config_file db.ini

### Run on one read pair

    nextflow run nextflow/variant_call.nf \
      -with-singularity clockwork_container.img \
      --ref_dir Reference.NC_000962.3/ \
      --reads_in1 reads_1.fastq --reads_in2 reads_2.fastq \
      --output_dir variant_call_out \
      --sample_name my_sample

# Clockwork 
## Pipelines for processing bacterial sequence data (Illumina only)  and variant calling

Note: these pipelines were developed for the [CRyPTIC](http://www.crypticproject.org/) project which studies _M. tuberculosis_, but in principle can be used on any bacteria.

Clockwork takes fastq input and outputs standard VCF files as output


Please see the [clockwork wiki page](https://github.com/iqbal-lab-org/clockwork/wiki) for documentation.

## Citation
If you use Clockwork for variant calling please cite the following paper, where it was introduced and benchmarked:  

"Minos: variant adjudication and joint genotyping of cohorts of bacterial genomes". Hunt et al. Genome Biol. 2022 Jul 5;23(1):147.
[doi: 10.1186/s13059-022-02714-x.](https://doi.org/10.1186/s13059-022-02714-x)

Whilst Minos is the novel part of the Clockwork variant calling pipeline, it uses several other tools. Please also cite:
* Minimap2 https://doi.org/10.1093/bioinformatics/bty191
* SAMtools https://doi.org/10.1093/gigascience/giab008
* Cortex https://doi.org/10.1038/ng.1028
* Trimmomatic: https://doi.org/10.1093/bioinformatics/btu170

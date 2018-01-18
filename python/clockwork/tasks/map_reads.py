from clockwork import read_map

def run(options):
    assert len(options.reads_files) % 2 == 0
    reads1_list = [filename for i, filename in enumerate(options.reads_files) if i % 2 == 0]
    reads2_list = [filename for i, filename in enumerate(options.reads_files) if i % 2 != 0]
    assert len(reads1_list) == len(reads2_list)
    read_map.map_reads_set(
        options.ref_fasta,
        reads1_list,
        reads2_list,
        options.outfile,
        rmdup=not options.unsorted_sam,
        read_group=('1', options.sample_name),
    )


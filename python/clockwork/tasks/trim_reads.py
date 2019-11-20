from clockwork import read_trim


def run(options):
    assert len(options.reads_files) % 2 == 0

    for i in range(0, len(options.reads_files), 2):
        out1 = options.outprefix + "." + str(int(i / 2)) + ".1.fq"
        out2 = options.outprefix + "." + str(int(i / 2)) + ".2.fq"
        read_trim.run_trimmomatic(
            options.reads_files[i],
            options.reads_files[i + 1],
            out1,
            out2,
            qual_trim=options.trimmo_qual,
        )

from clockwork import contam_remover

def run(options):
    remover = contam_remover.ContamRemover(
      options.metadata_tsv,
      options.bam_in,
      options.counts_out,
      options.reads_out_1,
      options.reads_out_2,
      contam_out_1=options.contam_out_1,
      contam_out_2=options.contam_out_2,
      no_match_out_1=options.no_match_out_1,
      no_match_out_2=options.no_match_out_2,
      done_file=options.done_file,
    )
    remover.run()


from clockwork import fake_contam_remover

def run(options):
    remover = fake_contam_remover.FakeContamRemover(
      options.reads_in_1,
      options.reads_in_2,
      options.reads_out_1,
      options.reads_out_2,
      options.counts_out,
    )
    remover.run()


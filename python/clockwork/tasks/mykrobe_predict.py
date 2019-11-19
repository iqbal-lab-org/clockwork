import sys
from clockwork import mykrobe


def run(options):
    panel = mykrobe.Panel(options.panel_dir)
    species = panel.metadata["species"]

    if panel.metadata["is_built_in"]:
        custom_probe_and_json = None
        panel_name = panel.metadata["name"]
    else:
        custom_probe_and_json = (panel.probes_fasta, panel.var_to_res_json)
        panel_name = None

    # This is for when we run the test in nextflow_mykrobe_predict_test.py.
    # Depending on sample, we want it to be resistant, susceptible,
    # or the run to fail
    if options.testing and options.sample_name.endswith("1_2"):
        sys.exit(1)

    mykrobe.run_predict(
        options.reads_files,
        options.output_dir,
        options.sample_name,
        species,
        panel=panel_name,
        custom_probe_and_json=custom_probe_and_json,
        unittest=options.testing,
    )

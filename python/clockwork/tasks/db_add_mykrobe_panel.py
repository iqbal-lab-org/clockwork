from clockwork import db


def run(options):
    database = db.Db(options.db_config_file)
    database.add_mykrobe_custom_panel(
        options.species,
        options.panel_name,
        options.reference_root,
        probes_fasta=options.probes_fasta,
        var_to_res_json=options.var_to_res_json,
    )

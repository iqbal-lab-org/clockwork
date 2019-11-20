from clockwork import db_maker


def run(options):
    dbm = db_maker.DbMaker(options.db_config_file)
    dbm.run()

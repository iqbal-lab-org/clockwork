import unittest
import configparser
import os
import pymysql
from clockwork import db_connection, db_maker, db_schema

modules_dir = os.path.dirname(os.path.abspath(db_maker.__file__))

class TestDbMaker(unittest.TestCase):
    def test_run(self):
        '''test run'''
        # use the test ini file. But we don't want the db line,
        # because we're making a new db
        ini_file = os.path.join(modules_dir, 'tests', 'data', 'db.ini')

        # in case database already exists
        try:
            db_connection.DbConnection(ini_file, destroy=True)
        except:
            pass
        # thorws error because database doesn't exist
        with self.assertRaises(db_connection.Error):
            db_connection.DbConnection(ini_file)


        dbm = db_maker.DbMaker(ini_file)
        dbm.run()

        # We'll just check that the database got created, and the
        # expected tables are in there. OTT to check complete
        # schema, methinks.
        dbc = db_connection.DbConnection(ini_file)
        cursor = dbc.connection.cursor()
        cursor.execute('USE ' + dbc.db)
        cursor.execute('show tables')
        got_tables = list(cursor.fetchall())
        got_tables.sort()
        expected_tables = [(x, ) for x in sorted(db_schema.tables)]
        self.assertEqual(expected_tables, got_tables)

        # check version got added
        got_rows = cursor.execute('SELECT * FROM Version;')
        expected_rows = [(db_schema.version)]
        dbc.close()
        db_connection.DbConnection(ini_file, destroy=True)

import unittest
import configparser
import os
import pymysql
from clockwork import db_connection

modules_dir = os.path.dirname(os.path.abspath(db_connection.__file__))
data_dir = os.path.join(modules_dir, "tests", "data", "db_connection")


class TestDbConnection(unittest.TestCase):
    def test_parse_config_file(self):
        """test _parse_config_file"""
        with self.assertRaises(Exception):
            db_connection.DbConnection._parse_config_file(
                os.path.join(data_dir, "config.missing_user.ini")
            )

        with self.assertRaises(Exception):
            db_connection.DbConnection._parse_config_file(
                os.path.join(data_dir, "config.no_db_login_header.ini")
            )

        got = db_connection.DbConnection._parse_config_file(
            os.path.join(data_dir, "config.good.ini")
        )
        expected = {
            "user": "username",
            "host": "hostname",
            "password": "pword",
            "db": "dbname",
            "port": None,
        }

        self.assertEqual(expected, got)

        got = db_connection.DbConnection._parse_config_file(
            os.path.join(data_dir, "config.good_with_port.ini")
        )
        expected = {
            "user": "username",
            "host": "hostname",
            "password": "pword",
            "db": "dbname",
            "port": 3306,
        }

        self.assertEqual(expected, got)

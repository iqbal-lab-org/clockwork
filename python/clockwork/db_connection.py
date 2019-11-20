import configparser
import pymysql


class Error(Exception):
    pass


class DbConnection:
    def __init__(self, ini_file, create=False, destroy=False, must_exist=False):
        login_info = DbConnection._parse_config_file(ini_file)
        db = None if create else login_info["db"]

        try:
            self.connection = pymysql.connect(
                host=login_info["host"],
                user=login_info["user"],
                password=login_info["password"],
                db=db,
                port=login_info["port"],
            )
        except:
            raise Error("Error connecting to database")

        if destroy:
            if must_exist:
                self.connection.cursor().execute("DROP DATABASE " + login_info["db"])
            else:
                self.connection.cursor().execute(
                    "DROP DATABASE IF EXISTS " + login_info["db"]
                )
            self.connection.commit()
            self.connection.close()
        elif create:
            self.connection.cursor().execute("CREATE DATABASE " + login_info["db"])
            self.connection.cursor().execute("USE " + login_info["db"])
            self.connection.commit()

        self.host = login_info["host"]
        self.user = login_info["user"]
        self.password = login_info["password"]
        self.db = login_info["db"]
        self.port = login_info["port"]

    @classmethod
    def _parse_config_file(cls, ini_file):
        config = configparser.ConfigParser()
        try:
            config.read(ini_file)
        except:
            raise Error("Error! [db_login] section not in file " + ini_file)

        expected_keys = {"user", "password", "host", "db"}
        got_keys = set(config["db_login"].keys())
        missing_keys = expected_keys.difference(got_keys)
        if len(missing_keys):
            raise Error(
                "Error! Missing these keys from [db_login] section: "
                + ",".join(list(missing_keys))
            )

        login_config = dict(config["db_login"])
        login_config["port"] = (
            int(login_config["port"]) if "port" in login_config else None
        )
        return login_config

    def commit(self):
        self.connection.commit()

    def close(self):
        self.connection.close()

    def reconnect(self):
        self.connection.commit()
        self.connection.close()
        self.connection = pymysql.connect(
            host=self.host, user=self.user, password=self.password, db=self.db
        )

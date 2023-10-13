from clockwork import db_connection, db_schema


class DbMaker:
    def __init__(self, ini_file):
        try:
            self.dbc = db_connection.DbConnection(ini_file, create=True)
        except:
            raise Exception("Error connecting to database")

    @classmethod
    def _create_table_command(cls, table, data_list, primary_key=None):
        command = ["CREATE TABLE", table, "("]

        for column, data_type in data_list:
            command.extend([column, data_type])

            if primary_key == column:
                command.append("AUTO_INCREMENT primary key")

            command.append(",")

        command.pop()
        command.append(")")
        return " ".join(command)

    def run(self):
        cursor = self.dbc.connection.cursor()
        for table_name, table_data in db_schema.tables.items():
            primary_key = db_schema.primary_keys.get(table_name, None)
            cursor.execute(
                DbMaker._create_table_command(
                    table_name, table_data, primary_key=primary_key
                )
            )
        cursor.execute("INSERT INTO Version VALUES (" + str(db_schema.version) + ")")
        self.dbc.commit()
        self.dbc.close()

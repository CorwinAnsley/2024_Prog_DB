import sqlite3
import os

CREATE_DB_SQLS = '/data/create_table_sqls/'

def create_db(db_filepath, sql_dir = CREATE_DB_SQLS):
    for sql in os.listdir(sql_dir):
        if os.isfile(os.join(sql_dir, sql)) and sql.endswith('.sql'):
            with open(os.join(sql_dir, sql)) as sql_file:
                stat = str(sql_file)
                print(stat)

if __name__ == '__main__':
    db_filepath = "test.db"
    db_connection = sqlite3.connection

    create_db(db_filepath)
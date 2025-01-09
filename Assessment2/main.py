import sqlite3
import os

from load_db import load_subjects

CREATE_DB_SQLS = './create_table_sqls/'

# Create the sqlite database if it does not exist and create the tables for all the entities
# using predefined sql scripts
def create_db(db_filepath, sql_dir = CREATE_DB_SQLS):
    # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    # Iterate through sql script
    for filepath in os.listdir(sql_dir):
        # Make this it
        if os.path.isfile(os.path.join(sql_dir, filepath)) and filepath.endswith('.sql'):
            with open(os.path.join(sql_dir, filepath)) as sql_file:
                sql = ''
                for line in sql_file:
                    sql += line#.strip()
                try:
                    db_cur.execute(sql)
                except:
                    pass
    
    db_connection.commit()
    db_cur.close()
    db_connection.close()


if __name__ == '__main__':
    db_filepath = "test.db"
    subjects_filepath = './data/Subject.csv'

    create_db(db_filepath)

    load_subjects(db_filepath, subjects_filepath)
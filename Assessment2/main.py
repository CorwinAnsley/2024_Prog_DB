import sqlite3
import os
import argparse
import matplotlib.pyplot as plt

import consts
from load_db import load_db

# Create the sqlite database if it does not exist and create the tables for all the entities
# using predefined sql scripts
def create_db(db_filepath, sql_dir = consts.CREATE_DB_SQLS):
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
                    sql += line
                try:
                    db_cur.execute(sql)
                except:
                    pass
    
    db_connection.commit()
    db_cur.close()
    db_connection.close()

# Takes database filepath and query number to execute one of the predefined queries
# Prints the results nicely and, for query 9, produces a plot
def query_db(db_filepath, query_no, sql_dir = consts.QUERY_DB_SQLS):
    # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    params = []
    sql = ''
    # Find the right query sql file
    for filepath in os.listdir(sql_dir):
        if filepath == f'query_{query_no}.sql':
            with open(os.path.join(sql_dir, filepath)) as sql_file:
                for line in sql_file:
                    sql += line
    

    db_cur.execute(sql, params)
    results = db_cur.fetchall()
    # Print out nicely
    print_q_result(results)

    # Create plot for query 9
    if query_no == 9:
        x = []
        y = []
        for row in results:
            x.append(row[0])
            y.append(row[1])
        
        plt.scatter(x, y)
        plt.xlabel('Age')
        plt.ylabel('BMI')
        plt.savefig('age_bmi_scatterplot.png')


# Takes the raw results and prints them out in a more readable fashion
def print_q_result(result):
    for row in result:
        s_row = []
        for i in range(len(row)):
            s_row.append(str(row[i]))
        row_str = '\t'.join(s_row)
        print(row_str)



    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='Database Coursework', description='Creates and manipulates a database.')
    parser.add_argument('--createdb',help='Create the database',action='store_true')
    parser.add_argument('--loaddb',help='load data into the database',action='store_true')
    parser.add_argument('--querydb',type=int,help='Specify a query to execute on the database',default=None)
    parser.add_argument('db_path', nargs='?', default=os.getcwd())

    args = parser.parse_args()
    db_filepath = str(args.db_path)
    createdb = args.createdb
    loaddb = args.loaddb
    querydb = args.querydb
    try:
        assert db_filepath.endswith('.db')
    except:
        raise Exception('Invalid database path')

    if createdb:
        create_db(db_filepath)

    if loaddb:
        load_db(db_filepath)
    
    if querydb:
        query_db(db_filepath, 9)
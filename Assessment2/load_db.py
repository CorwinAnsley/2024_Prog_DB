import sqlite3
import csv

import consts

def check_subject_entry(values):
    try:
        # Make sure there is an SubjectID
        values[0] = str(values[0])
        assert len(values[0]) > 0

        # Error check race value
        values[1] = str(values[1])
        try:
            assert len(values[1]) == 1
        except:
            values[1] = None

        # Error check sex value
        values[2] = str(values[2])
        try:
            assert values[2] in ['M', 'F']
        except:
            values[2] = None
        
        # Error check Age, BMI, SSPG
        for i in range(3,6):
            values[i] = float(values[i])
        
        # Error check IR-IS classification value
        values[6] = str(values[6])
        try:
            assert values[6] in ['IR', 'IS']
        except:
            values[6] = None
        
    except:
        raise Exception('Invalid entry')

def create_header_dict(headers, headers_order):
    header_dict = {}
    for i in range(len(headers)):
        for j in range(len(headers_order)):
            if headers[i] == headers_order[j]:
                header_dict[i] = j

def load_file_data(subjects_filepath, delimeter=',', headers_order = None):
    with open(subjects_filepath) as subject_file:
        reader = csv.reader(subject_file, delimiter=delimeter)
        headers = next(reader)
        if headers_order:
            header_dict = create_header_dict(headers, consts.SUBJECT_HEADERS_ORDER)
        
        file_data = []
        for row in reader:
            if headers_order:
                row_data = []
                for i in range(len(row)):
                    row_data[header_dict[i]] = row[i]
                file_data.append(row_data)
            else:
                file_data.append(row)
    
    return file_data
    
                
    # # create db connection
    # db_connection = sqlite3.connect(db_filepath)
    # db_cur = db_connection.cursor()



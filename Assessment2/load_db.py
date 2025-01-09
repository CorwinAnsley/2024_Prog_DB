import sqlite3
import csv

import consts

# Checks that an entry is valid to be inserted to the Subjct Table
def check_subject_entry(values, header_dict):
    try:
        #Ensure right amount of values
        assert len(values) == 7

        r_values = [None]*7
        # Make sure there is an SubjectID
        r_values[0] = str(values[header_dict[0]])
        assert len(r_values[0]) > 0

        # Error check race value
        r_values[1] = str(values[header_dict[1]])
        try:
            assert len(r_values[1]) == 1
        except:
            r_values[1] = None

        # Error check sex value
        r_values[2] = str(values[header_dict[2]])
        try:
            assert r_values[2] in ['M', 'F']
        except:
            r_values[2] = None
        
        # Error check Age, BMI, SSPG
        for i in range(3,6):
            try:
                r_values[i] = float(values[header_dict[i]])
            except:
                r_values[i] = None
        
        # Error check IR-IS classification value
        r_values[6] = str(values[header_dict[6]])
        try:
            assert r_values[6] in ['IR', 'IS']
        except:
            r_values[6] = None

        return r_values
    except Exception as err:
        # If any error occurs return False to be filtered out
        print(err)
        return False

# def check_subject_entry(values):
#     try:
#         #Ensure right amount of values
#         assert len(values) == 7

#         # Make sure there is an SubjectID
#         values[0] = str(values[0])
#         assert len(values[0]) > 0

#         # Error check race value
#         values[1] = str(values[1])
#         try:
#             assert len(values[1]) == 1
#         except:
#             values[1] = None

#         # Error check sex value
#         values[2] = str(values[2])
#         try:
#             assert values[2] in ['M', 'F']
#         except:
#             values[2] = None
        
#         # Error check Age, BMI, SSPG
#         for i in range(3,6):
#             values[i] = float(values[i])
        
#         # Error check IR-IS classification value
#         values[6] = str(values[6])
#         try:
#             assert values[6] in ['IR', 'IS']
#         except:
#             values[6] = None

#         return values
#     except:
#         return False

def check_visit_entry(values):
    try:
        #Ensure right amount of values
        assert len(values) == 2

        #Check SubjectID and VisitID
        values[0] = str(values[0])
        assert len(values[0]) > 0

        values[1] = str(values[1])
        assert len(values[1]) > 0 

        return values       
    except Exception as err:
        # If any error occurs return False to be filtered out
        print('aaaa')
        print(err)
        return False

def check_measurement_entry(values):
    try:
        #Ensure right amount of values
        assert len(values) == 4

        #Check Name, SubjectID and VisitID
        values[0] = str(values[0])
        assert len(values[0]) > 0

        values[1] = str(values[1])
        assert len(values[1]) > 0

        values[2] = str(values[2])
        assert len(values[2]) > 0

        try:
            values[3] = float(values[3])
        except:
            values[3] = None

        return values       
    except Exception as err:
        # If any error occurs return False to be filtered out
        print(err)
        return False

# Create dict with the key of table entry and value of position in the rows of headers
def create_header_dict(headers, headers_order):
    header_dict = {}
    for i in range(len(headers)):
        for j in range(len(headers_order)):
            if headers[i] == headers_order[j]:
                header_dict[j] = i
    return header_dict

def load_file_data(subjects_filepath, delimeter=',', headers_order = None):
    with open(subjects_filepath, 'r', encoding='utf-8-sig') as subject_file:
        reader = csv.reader(subject_file, delimiter=delimeter)
        headers = next(reader)
        if headers_order:
            header_dict = create_header_dict(headers, consts.SUBJECT_HEADERS_ORDER)
        
        file_data = []
        row_len = len(header_dict)
        for row in reader:
            if headers_order:
                row_data = [None] * row_len
                for i in range(len(row)):
                    try:
                        row_data[header_dict[i]] = row[i]
                    except Exception as err:
                        print(err)
                file_data.append(row_data)
            else:
                file_data.append(row)
    
    return file_data

# Insert entries into table, ignoring duplicates
def insert_entries(data, table, db_cur, num_values):
    sql = f'INSERT OR IGNORE INTO {table} VALUES('
    for i in range(num_values-1):
        sql += '?, '
    
    sql += '?);'
    #print(next(data))
    db_cur.executemany(sql, data)


# Load all entries from a file to table in database
def load_entries(db_filepath, table_name, entries_filepath, headers_order, delimeter=','):
    # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    #file_data = load_file_data(subjects_filepath, delimeter=',',headers_order=consts.SUBJECT_HEADERS_ORDER)

    # Open file with entries
    with open(entries_filepath, 'r', encoding='utf-8-sig') as entries_file:
        reader = csv.reader(entries_file, delimiter=delimeter)
        # Get headers and create dict to shift columns as expected for db table
        headers = next(reader)
        header_dict = create_header_dict(headers, headers_order)

        # Each entry is checked for errors and corrected, or filtered out
        entries_data = filter(bool, map(lambda entry : check_subject_entry(entry, header_dict), reader))
        insert_entries(entries_data, table_name, db_cur, len(headers_order))

    db_connection.commit()
    db_cur.close()
    db_connection.close()


def load_measurements_visits(db_filepath, tablename_measures, tablename_visits, tablename_measure_names, type, entries_filepath, delimeter='\t'):
        # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    # Open file with entries
    with open(entries_filepath, 'r', encoding='utf-8-sig') as entries_file:
        reader = csv.reader(entries_file, delimiter=delimeter)
        # Get headers, check they are all valid and use to populate MeasurementName table
        headers = next(reader)
        measurement_names = []
        for header in headers:
            try:
                header = str(header)
                if len(header) > 0:
                    measurement_names.append([header,type])
            except:
                pass
        
        insert_entries(measurement_names, tablename_measure_names, db_cur, 2)

        #parse_measurement_row(db_cur, tablename_measures, tablename_visits, next(reader), headers)
        #map(lambda row : parse_measurement_row(db_cur, tablename_measures, tablename_visits, row, headers), reader)
        for row in reader:
            parse_measurement_row(db_cur, tablename_measures, tablename_visits, row, headers)
    
    db_connection.commit()
    db_cur.close()
    db_connection.close()

def parse_measurement_row(db_cur, tablename_measures, tablename_visits, row, headers):
    try:
        SubjectID, VisitID = row[0].split('-')
        assert check_visit_entry([VisitID, SubjectID])
        insert_entries([[VisitID,SubjectID]], tablename_visits, db_cur, 2)
        row_data = []
        for i in range(1,len(row)):
            values = [headers[i],VisitID,SubjectID,row[i]]
            values = check_measurement_entry(values)
            if values:
                row_data.append([headers[i],VisitID,SubjectID,row[i]])
        insert_entries(row_data, tablename_measures, db_cur, 4)
    except Exception as err:
        print('n')
        print(err)
        return
        

def load_db(db_filepath):
    # Load Suject entries
    #load_entries(db_filepath, 'Subject', consts.SUBJECT_FILEPATH, consts.SUBJECT_HEADERS_ORDER, delimeter=',')

    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Transcriptome', consts.TRANSCRIPTOME_MEASURES_FILEPATH)
    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Proteome', consts.PROTEOME_MEASURES_FILEPATH)
    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Metabolome', consts.METABOLOME_MEASURES_FILEPATH)
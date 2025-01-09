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
        print(err)
        return False

# Error check measurements
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

# Error check peaks and metabolomes
def check_peak_metab_entry(values, header_dict):
    try:
        print(header_dict)
        r_values = [None]*7

        #Check SubjectID and VisitID
        r_values[0] = str(values[header_dict[0]])
        assert len(r_values[0]) > 0

        for i in range(1,5):
            r_values[i] = str(values[header_dict[i]])
        
        r_values[5] = str(values[header_dict[5]])
        assert len(r_values[1]) > 0

        print('a')
        #Copy Metabolite nname for PeakID table
        r_values[6] = r_values[0]

        return r_values
    except Exception as err:
        # If any error occurs return False to be filtered out
        print(err)
        print('e')
        return False

# Create dict with the key of table entry and value of position in the rows of headers
def create_header_dict(headers, headers_order):
    header_dict = {}
    for i in range(len(headers)):
        for j in range(len(headers_order)):
            if headers[i] == headers_order[j]:
                header_dict[j] = i
    return header_dict

# Insert entries into table, ignoring duplicates
def insert_entries(data, table, db_cur, num_values):
    sql = f'INSERT OR IGNORE INTO {table} VALUES('
    for i in range(num_values-1):
        sql += '?, '
    
    sql += '?);'

    db_cur.executemany(sql, data)

# Special function to insert both Peaks and Metabolites
# Two sql statements are made for each then the data is inserted from each row
def insert_peak_metab_entries(data, table, table2, db_cur, num_values, table2_num_values):
    table1_num_values = num_values - table2_num_values
    sql1 = f'INSERT OR IGNORE INTO {table} VALUES('
    for i in range(table1_num_values-1):
        sql1 += '?, '
    
    sql1 += '?);'

    sql2 = f'INSERT OR IGNORE INTO {table2} VALUES('
    for i in range(table1_num_values,num_values-1):
        sql2 += '?, '
    
    sql2 += '?);'
    #sql = f'{sql1}\n{sql2}'
  
    for row in data:
        print(row[:table1_num_values])
        db_cur.execute(sql1,row[:table1_num_values])
        db_cur.execute(sql2,row[table1_num_values:])
    #db_cur.executemany(sql, data)


# Load all entries from a file to table in database
# Used for Subjects file as well as Peak/Metabolite file
def load_entries(db_filepath, table_name, entries_filepath, headers_order, check_func, delimeter=',', table2_name = '', table2_num_values = 0):
    # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    # Open file with entries
    with open(entries_filepath, 'r', encoding='utf-8-sig') as entries_file:
        reader = csv.reader(entries_file, delimiter=delimeter)
        # Get headers and create dict to shift columns as expected for db table
        headers = next(reader)
        header_dict = create_header_dict(headers, headers_order)

        # Each entry is checked for errors and corrected, or filtered out
        entries_data = filter(bool, map(lambda entry : check_func(entry, header_dict), reader))
        # For Peaks/Metabolites
        if table2_name != '':
            insert_peak_metab_entries(entries_data, table_name, table2_name, db_cur, len(headers_order)+1, table2_num_values=table2_num_values)
        # For Subjects
        else:
            insert_entries(entries_data, table_name, db_cur, len(headers_order))

    db_connection.commit()
    db_cur.close()
    db_connection.close()

# Load Measurements and visits from file
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

        # Iterate through each row and load into the appropiate tables
        for row in reader:
            parse_measurement_row(db_cur, tablename_measures, tablename_visits, row, headers)
    
    # Commit changes
    db_connection.commit()
    db_cur.close()
    db_connection.close()

# Parse each row in measurement file for data to insert
def parse_measurement_row(db_cur, tablename_measures, tablename_visits, row, headers):
    try:
        # Obtain SubjectID and VisitID by splitting first item
        SubjectID, VisitID = row[0].split('-')

        # Check the visit is valid
        assert check_visit_entry([VisitID, SubjectID])

        #Insert the visit
        insert_entries([[VisitID,SubjectID]], tablename_visits, db_cur, 2)

        # Iterate through each item, check it is valid measurement and add to row_data
        row_data = []
        for i in range(1,len(row)):
            values = [headers[i],VisitID,SubjectID,row[i]]
            values = check_measurement_entry(values)
            if values:
                row_data.append([headers[i],VisitID,SubjectID,row[i]])
        
        #Insert all valid measurements
        insert_entries(row_data, tablename_measures, db_cur, 4)
    except:
        return

def load_db(db_filepath):
    # Load Subject entries
    load_entries(db_filepath, 'Subject', consts.SUBJECT_FILEPATH, consts.SUBJECT_HEADERS_ORDER, check_subject_entry, delimeter=',')

    # Load Measurements, MeasurementNames and Visits
    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Transcriptome', consts.TRANSCRIPTOME_MEASURES_FILEPATH)
    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Proteome', consts.PROTEOME_MEASURES_FILEPATH)
    load_measurements_visits(db_filepath, 'Measurement', 'Visit', 'MeasurementName', 'Metabolome', consts.METABOLOME_MEASURES_FILEPATH)

    # Load Peaks and Metabolites
    load_entries(db_filepath, 'Metabolite', consts.PEAKS_METABOLITES_FILEPATH, consts.PEAKS_METABOLITES_HEADERS_ORDER,check_peak_metab_entry, delimeter=',',table2_name='Peak', table2_num_values=2)
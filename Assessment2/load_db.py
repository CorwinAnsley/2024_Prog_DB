import sqlite3
import csv

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

def load_subjects(db_filepath, subjects_filepath):
    # create db connection
    db_connection = sqlite3.connect(db_filepath)
    db_cur = db_connection.cursor()

    with open(subjects_filepath) as subject_file:
        spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')

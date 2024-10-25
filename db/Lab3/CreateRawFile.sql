
CREATE TABLE RawFile (
  SampleID VARCHAR(255) NOT NULL,
  RawFileID VARCHAR(255) NOT NULL,
  OmicsType VARCHAR(20),
  Filename VARCHAR(255),
  Timestamp INT,
  BatchNum INT,
  RunType VARCHAR(255),
  
  PRIMARY KEY (RawFileID)
  FOREIGN KEY (SampleID) REFERENCES Patient(SampleID)
);
-- Create a table that corresponds to the Patient entity
CREATE TABLE Patient (
  SampleID VARCHAR(255) NOT NULL,
  SampleName VARCHAR(255),
  Age INT,
  ICU INT,
  Gender VARCHAR(255),
  IsCovid INT,
  OnVentilator INT,
  HFD_45 INT,
  WHO_Ordinal INT,
  APACHE_II INT,
  SOFA INT,
  CharlsonScore INT,
  
  PRIMARY KEY (SampleID)
);
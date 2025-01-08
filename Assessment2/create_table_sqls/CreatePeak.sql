-- Create table for the Peak entity

 CREATE TABLE Peak(
  PeakID VARCHAR(255) NOT NULL,
  MetabolomeName VARCHAR(255) NOT NULL,
  FOREIGN KEY (MetabolomeName) REFERENCES Metabolome(Name),
  PRIMARY KEY (MetabolomeName, PeakID)
); 
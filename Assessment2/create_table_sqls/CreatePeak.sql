-- Create table for the Peak entity

 CREATE TABLE Peak(
  PeakID VARCHAR(255) NOT NULL,
  MetaboliteName VARCHAR(255) NOT NULL,
  FOREIGN KEY (MetaboliteName) REFERENCES Metabolite(Name),
  PRIMARY KEY (MetaboliteName, PeakID)
); 
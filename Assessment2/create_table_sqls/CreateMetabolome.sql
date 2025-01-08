-- Create table for the Metabolome entity

 CREATE TABLE Metabolome(
  Name VARCHAR(255) NOT NULL,
  Pathway VARCHAR(255),
  Chemical_Class VARCHAR(255),
  HMDB VARCHAR(255),
  KEGG VARCHAR(255),
  PRIMARY KEY (Name)
); 
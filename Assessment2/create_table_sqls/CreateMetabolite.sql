-- Create table for the Metabolite entity

 CREATE TABLE Metabolite(
  Name VARCHAR(255) NOT NULL,
  Pathway VARCHAR(255),
  Chemical_Class VARCHAR(255),
  HMDB VARCHAR(255),
  KEGG VARCHAR(255),
  PRIMARY KEY (Name)
); 
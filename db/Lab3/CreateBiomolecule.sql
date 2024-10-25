CREATE TABLE Biomolecule(
  BiomoleculeID VARCHAR(255) NOT NULL,
  StandardisedName VARCHAR(255),
  OmicsType VARCHAR(20),
  KeepNote VARCHAR(255),
  
  PRIMARY KEY (BioMoleculeID)
);

CREATE TABLE TestResult(
  TestName VARCHAR(255) NOT NULL,
  BiomoleculeID VARCHAR(255) NOT NULL,
  pValue DECIMAL(10, 5),
  qValue DECIMAL(10, 5),
  
  PRIMARY KEY (BioMoleculeID, TestName)
  FOREIGN KEY (BioMoleculeID) REFERENCES BioMolecule(BioMoleculeID)
);
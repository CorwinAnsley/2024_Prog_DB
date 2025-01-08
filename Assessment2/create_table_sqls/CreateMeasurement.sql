-- Create table for the Measurement entity

 CREATE TABLE Measurement(
  Name VARCHAR(255) NOT NULL,
  VisitID VARCHAR(255) NOT NULL,
  SubjectID VARCHAR(255) NOT NULL,
  Value DECIMAL (13,6),
  FOREIGN KEY (VisitID) REFERENCES Visit(VisitID),
  FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID),
  PRIMARY KEY (SubjectID, VisitID, Name)
);  
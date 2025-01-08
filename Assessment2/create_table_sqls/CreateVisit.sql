-- Create table for the Visit entity

 CREATE TABLE Visit(
  VisitID VARCHAR(255) NOT NULL,
  SubjectID VARCHAR(255) NOT NULL,
  FOREIGN KEY (SubjectID) REFERENCES Subject(SubjectID),
  PRIMARY KEY (SubjectID, VisitID)
);  
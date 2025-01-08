-- Create a table for the Subject entity 

 CREATE TABLE Subject(
  SubjectID VARCHAR(255) NOT NULL,
  Race CHAR(1),
  Sex CHAR(1), -- 'M' or 'F'
  Age DECIMAL (5,2),
  BMI DECIMAL (5,2),
  SSPG DECIMAL (6,2),
  IR_IS_classification CHAR(2), -- 'IR' or 'IS'
  PRIMARY KEY (SubjectID)
);  
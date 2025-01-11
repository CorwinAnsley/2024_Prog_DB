SELECT MAX(Measurement.Value)
FROM Measurement
INNER JOIN Subject
	ON Measurement.SubjectID = Subject.SubjectID
WHERE Measurement.Name = 'A1BG' AND Subject.SubjectID = 'ZOZOW1T' 
 
SELECT Subject.SubjectID
FROM Subject
INNER JOIN Measurement
	ON Subject.SubjectID = Measurement.SubjectID
INNER JOIN MeasurementName
	ON Measurement.Name = MeasurementName.Name
WHERE Subject.IR_IS_classification = 'IS' AND MeasurementName.Type = 'Metabolome'
GROUP BY Subject.SubjectID

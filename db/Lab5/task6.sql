SELECT FileName
FROM RawFile
	INNER JOIN Patient
		ON RawFile.SampleID = Patient.SampleID
	WHERE RawFile.OmicsType = 'metabolomics' AND
		Patient.SampleName = 'COVID_01'
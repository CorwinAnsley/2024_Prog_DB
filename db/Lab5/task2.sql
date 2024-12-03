SELECT SampleName
FROM Patient
WHERE EXISTS (
	SELECT Patient.SampleName
	FROM Patient
		INNER JOIN RawFile
		ON Patient.SampleID = RawFile.SampleID
	WHERE RawFile.OmicsType = 'transcriptomics'	
);
	
SELECT Patient.SampleName, Rawfile.FileName, RawFile.OmicsType
FROM Patient
	INNER JOIN RawFile
		ON Patient.SampleID = RawFile.SampleID 
WHERE Patient.SampleName IN ('COVID_02', 'COVID_04', 'COVID_06', 'COVID_8','COVID_10')
ORDER BY Patient.SampleName, RawFile.FileName;
SELECT COUNT(SampleID), Gender
FROM Patient
WHERE CovidStatus = '1'
GROUP BY Gender

SELECT COUNT(SampleID), Gender, CovidStatus
FROM Patient
GROUP BY Gender, CovidStatus
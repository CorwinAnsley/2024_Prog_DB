SELECT COUNT(SampleID), AVG(AGE), Gender, CovidStatus
FROM Patient
GROUP BY Gender, CovidStatus
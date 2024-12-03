SELECT AVG(Patient.HDF45) AS Mean_HDF45,
	MIN(Patient.HDF45) AS Lowest_HDF45,
	Max(Patient.HDF45) AS Highest_HDF45
FROM Patient
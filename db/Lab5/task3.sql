SELECT *
FROM Patient
WHERE Patient.CovidStatus = '0' AND
	Patient.Age > (SELECT AVG(Age)
		FROM Patient
		WHERE Patient.CovidStatus = '1'
	)
ORDER BY Age DESC;
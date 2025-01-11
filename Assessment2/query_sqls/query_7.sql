SELECT Pathway, COUNT(Pathway)
FROM Metabolite
INNER JOIN Peak	
 ON Metabolite.Name = Peak.MetaboliteName
GROUP BY Pathway
HAVING COUNT(Pathway) > 9
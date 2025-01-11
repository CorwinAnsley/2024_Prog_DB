SELECT Metabolite.KEGG
FROM Metabolite
INNER JOIN Peak
	ON Peak.MetaboliteName = Metabolite.Name
WHERE Peak.PeakID IN ('nHILIC_121.0505_3.5', 'nHILIC_130.0872_6.3', 'nHILIC_133.0506_2.3','nHILIC_133.0506_4.4')
GROUP BY Metabolite.KEGG
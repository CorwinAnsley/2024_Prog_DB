SELECT BioMolecule.StandardisedName
		FROM BioMolecule
		INNER JOIN BiomoleculeMeasurement
			ON BioMolecule.BioMoleculeID = BiomoleculeMeasurement.BioMoleculeID
		WHERE BiomoleculeMeasurement.NormalisedAbundance = 0
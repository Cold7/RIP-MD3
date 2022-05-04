import MDAnalysis as mda
def save_selection(selection, pdbName, outName):
	try:
		universe = mda.Universe(pdbName)
		#selecting atoms to work
		sel = universe.select_atoms(selection)
		writerPDB=mda.Writer(outName)
		writerPDB.write(sel)	
		return outName
	except:
		pass
	return False

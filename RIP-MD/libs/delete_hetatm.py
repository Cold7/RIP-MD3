import multiprocessing as mp
import MDAnalysis



def parallelDelete(dictFramePDB):
	pdb=dictFramePDB[1]
	pdbName=(pdb.split("/"))[-1]

	universe=MDAnalysis.Universe(pdb)
	selection=universe.select_atoms("protein") #in v<0.11 we need to use select_atoms. MODIFIED.
	
	writerPDB=MDAnalysis.Writer(Output+"/.temp_RIP-MD/"+str(pdbName[0:-4])+"_nohetatm.pdb")
	writerPDB.write(selection)
#	os.system("sed -i \"s/HSD\|HSP\|HSE/HIS/g\" "+Output+"/.temp_RIP-MD/"+str(pdbName[0:-4])+"_nohetatm.pdb") #we rename HSD or HSP to HIS because DSSP does not recognice those amino acids

	return (dictFramePDB[0],Output+"/.temp_RIP-MD/"+pdbName[0:-4]+"_nohetatm.pdb")
	
########################################################################
##
## function that create a file based on a pdb and delete all hetatm
##
########################################################################

def delete(pdbDict,output,nproc):
	global Output
	Output=output
	pool=mp.Pool(processes=int(nproc)) #for multiprocessing	
	listFrames=list(pool.map(parallelDelete,pdbDict.items(),chunksize=1))
	pool.close()
	DictPDB={}
	for lista in listFrames:
		DictPDB[lista[0]]=lista[1]
	listFrames=None
	return DictPDB

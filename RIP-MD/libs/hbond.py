import multiprocessing as mp # to parallelize jobs
import sys
import MDAnalysis
import MDAnalysis.analysis.hydrogenbonds as hbonds
import warnings
warnings.filterwarnings("ignore")


def parallelParser(pdbFrameName): #arg: list with universe, start and stop frame, distance and angle to calculate hbond
	sys.stdout = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')				
	sys.stderr = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')
	
	edges_file=open(outputFolder+"/RIP-MD_Results/Edges/"+pdbFrameName[0]+".edges","a") #0 contains the frame_number
	edges_file.write("Hydrogen Bonds\nSource Node\tTarget Node\tEdge Name\tDistance\tAngle\n")
	u = MDAnalysis.Universe(pdbFrameName[1])
	h = hbonds.HydrogenBondAnalysis(u,donors_sel="protein",hydrogens_sel="protein",acceptors_sel="protein", d_a_cutoff=Distance, d_h_a_angle_cutoff=Angle)
	h.run()


	# OLD VERSION
	# h.generate_table()
	#
	# for hbond in h.table:
	#
	# 	node1=str(hbond[5])+":"+dictAtom[hbond[1]]+":"+str(hbond[6])
	# 	node2=str(hbond[8])+":"+dictAtom[hbond[2]]+":"+str(hbond[9])
	# 	if (node1 in nodes) and (node2 in nodes):
	# 		name=node1+":"+str(hbond[7])+"-(HB)-"+node2+":"+str(hbond[10])
	# 		distance=str(round(float(str(hbond[11])),3))
	# 		angle=str(round(float(str(hbond[12])),3))
	# 		edges_file.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+name+"\t"+distance+"\t"+angle+"\n")
	# edges_file.write("\n")
	# edges_file.close()

	#NEW VERSION 1
	# for hbond in h.results.hbonds:
	# 	donor = u.select_atoms("index " + str(int(hbond[1])))
	# 	acceptor = u.select_atoms("index " + str(int(hbond[3])))
	#
	# 	node1 = str(donor[0].resname) + ":" + dictAtom[int(hbond[1])+1] + ":" + str(donor[0].resid)
	# 	node2 = str(acceptor[0].resname) + ":" + dictAtom[int(hbond[3])+1] + ":" + str(acceptor[0].resid)
	# 	if (node1 in nodes) and (node2 in nodes):
	# 		name=node1+":"+str(donor[0].type)+"-(HB)-"+node2+":"+str(acceptor[0].type)
	# 		distance=str(round(float(str(hbond[4])),3))
	# 		angle=str(round(float(str(hbond[5])),3))
	# 		edges_file.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+name+"\t"+distance+"\t"+angle+"\n")
	#
	# edges_file.write("\n")
	# edges_file.close()


	#NEW VERSION 2
	for hbond in h.results.hbonds:
		donor = u.atoms[int(hbond[1])]
		acceptor = u.atoms[int(hbond[3])]

		node1 = str(donor.resname) + ":" + str(donor.segid) + ":" + str(donor.resid)
		node2 = str(acceptor.resname) + ":" + str(acceptor.segid) + ":" + str(acceptor.resid)
		if (node1 in nodes) and (node2 in nodes):
			name=node1+":"+str(donor.type)+"-(HB)-"+node2+":"+str(acceptor.type)
			distance=str(round(float(str(hbond[4])),3))
			angle=str(round(float(str(hbond[5])),3))
			edges_file.write(str(nodes[node1])+"\t"+str(nodes[node2])+"\t"+name+"\t"+distance+"\t"+angle+"\n")

	edges_file.write("\n")
	edges_file.close()



#//////////////////////////////////////////////
#/
#/ main function that will call a parallel
#/ function. this parallel function will
#/ parse and return nodes and attrs for
#/ Hbonds
#/
#/////////////////////////////////////////////
def parse(pdbDict, OutputFolder, nproc, Nodes, RIPMDPath, Dist, Ang):
	global nodes, outputFolder, ripmdPath, Distance, Angle, dictAtom
	nodes=Nodes
	outputFolder=OutputFolder
	ripmdPath=RIPMDPath
	Distance=float(Dist)
	Angle=float(Ang)
	sys.stdout = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')	
	sys.stderr = open(outputFolder+"/.temp_RIP-MD/temp_hbond", 'a')

	# u = MDAnalysis.Universe(pdbDict["frame_0"])
	# dictAtom={}
	# for atom in u.atoms:
	# 	dictAtom[atom.id]=atom.segid
	chains=[]

	pool=mp.Pool(processes=int(nproc)) #for multiprocessing
	listOfAttr=list(pool.map(parallelParser,pdbDict.items(),chunksize=1))
	pool.close()
	sys.stdout = open(outputFolder+"/output.log", 'a')
	sys.stderr = open(outputFolder+"/output.log", 'a')

	return

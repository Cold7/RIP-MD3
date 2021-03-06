import networkx as nx
import saveGraph as sg
import multiprocessing as mp

def parallelGraph(arguments):
	pattern, toCompute, graphName=arguments
	graph=nx.MultiGraph()
	graph.add_nodes_from(nodeGraph.nodes(data=True))
	if toCompute:
		edges=open(folder+"/RIP-MD_Results/Edges/frame_0.edges","r")
		for line in edges:
			if pattern in line:
				splittedLine=line[:-1].split("\t")
				if pattern=="(Ca)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Alpha Carbon Distance",distance=str(splittedLine[3]))
				elif pattern=="(HB)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Hydrogen Bond",distance=str(splittedLine[3]), angle=str(splittedLine[4]))						
				elif pattern=="(salt)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Salt Bridge",distance=str(splittedLine[3]))
				elif pattern=="(SS)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Disulfide Bond",distance=str(splittedLine[3]), Angle=str(splittedLine[4]))											
				elif pattern=="(cation-pi)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Cation-Pi",distance=str(splittedLine[3]), orientation=str(splittedLine[4]))
				elif pattern=="(pi-pi)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Pi-Pi",distance=str(splittedLine[3]), orientation=str(splittedLine[4]))
				elif pattern=="(Arg-Arg)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Arg-Arg",distance=str(splittedLine[3]))													
				elif pattern=="(coulomb)":
					edgeName=(((splittedLine[2].replace("[","")).replace("]","")).replace("'","")).replace(" ","")
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=edgeName,interaction="Coulomb",distance=str(splittedLine[3]),potential=str(splittedLine[4])+" kJ/mole")						
				elif pattern=="(vdW)":
					graph.add_edge(nodeDict[splittedLine[0]],nodeDict[splittedLine[1]], name=splittedLine[2],interaction="Van der Waals",distance=str(splittedLine[3]), potential=str(splittedLine[4])+" kJ/mole")	
		edges.close()

	# removing nodes without adjacent edges
	listToRemove = []
	for node in graph.nodes():
		if graph.degree[node] == 0:
			listToRemove.append(node)
	graph.remove_nodes_from(listToRemove)

	if len(graph.nodes())>1:
			
		#now we will save the graph in a specific format
		if gFormat.lower()=="gml":
			sg.gml(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".gml")
		if gFormat.lower()=="edgelist":
			sg.edgelist(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".edgelist")
		if gFormat.lower()=="gefx":
			sg.gexf(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".gexf")	
		if gFormat.lower()=="graphmL":
			sg.graphml(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".graphml")
		if gFormat.lower()=="pajek":
			sg.pajek(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".pajek")	
	return graph

def MDProperties(parameters):
	ID,node=parameters
	sasa=0
	cont=0
	eightClass={".":0.0,"S":0.0,"T":0.0,"C":0.0,"H":0.0,"G":0.0,"I":0.0,"E":0.0,"B":0.0}
	threeClass={"C":0.0,"H":0.0,"E":0.0}
	for item in pdbNames.items():
		attrFile=open(folder+"/RIP-MD_Results/Attrs/"+item[0]+".attr","r")

		for line in attrFile:
			if line[:-1]!="Node ID\tSASA\t8-class sec. struct.\t3-class sec. stuct." and line!="\n":
				aux=line[:-1].split("\t") #pos. 0=id, 1=sasa, 2= sec struct 8 class, 3: sec struct 3 class
				if aux[0]==ID:
					sasa+=int(aux[1])
					eightClass[aux[2]]+=1
					threeClass[aux[3]]+=1
		attrFile.close()			
		cont+=1
	for Class in eightClass.items():
		eightClass[Class[0]]=(Class[1]*100)/cont
	for Class in threeClass.items():
		threeClass[Class[0]]=(Class[1]*100)/cont		
	sasa/=cont
	return [node,eightClass,threeClass,sasa]	

def MDGraph(arguments):
	
	pattern, toCompute, graphName=arguments
	graph=nx.MultiGraph()
	graph.add_nodes_from(nodeGraph.nodes(data=True))
	interactions={}
	if toCompute:
		for frame in pdbNames.items():
			edges=open(folder+"/RIP-MD_Results/Edges/"+frame[0]+".edges","r")
			for line in edges:
				if pattern in line:
					splittedLine=line[:-1].split("\t")
					####################################################
					##
					## for Calpha, salt and args
					##
					####################################################
					if pattern=="(Ca)" or pattern=="(salt)" or pattern=="(Arg-Arg)":
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							interactions[splittedLine[2]][2]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),1.0] #name, dist and time
					####################################################
					##
					## for H bonds
					##
					####################################################
					if pattern=="(HB)":
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							interactions[splittedLine[2]][2]+=float(splittedLine[4])
							interactions[splittedLine[2]][3]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),float(splittedLine[4]),1.0] #name, dist, angle and time
					####################################################
					##
					## for disulfide bonds
					##
					####################################################
					if pattern=="(SS)":
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							interactions[splittedLine[2]][2]+=float(splittedLine[4])
							interactions[splittedLine[2]][3]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),float(splittedLine[4]),1.0] #name, dist, diedralAngle and time
					####################################################
					##
					## for cation pi or pi-pi
					##
					####################################################
					if pattern=="(cation-pi)" or pattern=="(pi-pi)" : 
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							#getting the orientation dict
							orientationDict=interactions[splittedLine[2]][2]
							if splittedLine[4] in orientationDict:
								interactions[splittedLine[2]][2][splittedLine[4]]+=1
							else:
								interactions[splittedLine[2]][2][splittedLine[4]]=1
							interactions[splittedLine[2]][3]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),{splittedLine[4]:1.0},1.0] #name, dist, {orientation:number} and time
                
					####################################################
					##
					## for coulomb
					##
					####################################################
					if pattern=="(coulomb)":
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							interactions[splittedLine[2]][2]+=float(splittedLine[4])
							interactions[splittedLine[2]][3]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),float(splittedLine[4]),1.0] #name, dist, potential and time								
				
					####################################################
					##
					## for VDW
					##
					####################################################
					if pattern=="(vdW)":
						if splittedLine[2] in interactions:
							interactions[splittedLine[2]][1]+=float(splittedLine[3])
							interactions[splittedLine[2]][2]+=float(splittedLine[4])
							interactions[splittedLine[2]][3]+=1.0
						else:
							interactions[splittedLine[2]]=[splittedLine[2],float(splittedLine[3]),float(splittedLine[4]),1.0] #name, dist, potential and time							

			edges.close()

		

	if pattern=="(Ca)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][2]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(Ca)-")
				graph.add_edge(splitted[0],splitted[1], name=key,interaction="Alpha Carbon Distance",avgDistance=str(round(((interactions[key][1])/interactions[key][2]),3)), time=str(round(((interactions[key][2]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(salt)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][2]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(salt)-")
				node1=splitted[0].split(":")
				node1=node1[0]+":"+node1[1]+":"+node1[2]
				node2=splitted[1].split(":")
				node2=node2[0]+":"+node2[1]+":"+node2[2]
				graph.add_edge(node1,node2, name=key,interaction="Salt Bridge",avgDistance=str(round(((interactions[key][1])/interactions[key][2]),3)), time=str(round(((interactions[key][2]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(Arg-Arg)":	
		keys=list(interactions.keys())
		for key in keys:
			if (interactions[key][2]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(Arg-Arg)-")
				graph.add_edge(splitted[0],splitted[1], name=key,interaction="Arg-Arg",avgDistance=str(round(((interactions[key][1])/interactions[key][2]),3)), time=str(round(((interactions[key][2]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(HB)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(HB)-")
				node1=splitted[0].split(":")
				node1=node1[0]+":"+node1[1]+":"+node1[2]
				node2=splitted[1].split(":")
				node2=node2[0]+":"+node2[1]+":"+node2[2]
				graph.add_edge(node1,node2, name=key,interaction="Hydrogen Bond",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), avgAngle=str(round(((interactions[key][2])/interactions[key][3]),3)), time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(SS)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(SS)-")
				graph.add_edge(splitted[0],splitted[1], name=key,interaction="Disulfide Bridge",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), avgAngle=str(round(((interactions[key][2])/interactions[key][3]),3)), time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(cation-pi)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(cation-pi)-")
				orientations=interactions[key][2]
				total=0
				for item in orientations.items():
					total+=item[1]
				for item in orientations.items():
					orientations[item[0]]=str(round((orientations[item[0]]*100)/total,3))+"%"
				graph.add_edge(splitted[0],splitted[1], name=key,interaction="Cation-Pi",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), orientation=str(orientations)[1:-1], time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
	if pattern=="(pi-pi)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(pi-pi)-")
				orientations=interactions[key][2]
				total=0
				for item in orientations.items():
					total+=item[1]
				for item in orientations.items():
					orientations[item[0]]=str(round((orientations[item[0]]*100)/total,3))+"%"
				graph.add_edge(splitted[0],splitted[1], name=key,interaction="Pi-Pi",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), orientation=str(orientations)[1:-1], time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]
				
				
	if pattern=="(coulomb)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(coulomb)-")
				n1=splitted[0].split(":")
				n2=splitted[1].split(":")
				nod1=n1[0]+":"+n1[1]+":"+n1[2]
				nod2=n2[0]+":"+n2[1]+":"+n2[2]
				graph.add_edge(nod1,nod2, name=(((key.replace("[","")).replace("]","")).replace("'","")).replace(" ",""),interaction="Coulomb",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), avgPotential=str(round(((interactions[key][2])/interactions[key][3]),3))+" kJ/mole", time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]

	if pattern=="(vdW)":
		keys = list(interactions.keys())
		for key in keys:
			if (interactions[key][3]*100)/len(pdbNames)>=float(time):
				splitted=key.split("-(vdW)-")
				n1=splitted[0].split(":")
				n2=splitted[1].split(":")
				nod1=n1[0]+":"+n1[1]+":"+n1[2]
				nod2=n2[0]+":"+n2[1]+":"+n2[2]
				graph.add_edge(nod1,nod2, name=(((key.replace("[","")).replace("]","")).replace("'","")).replace(" ",""),interaction="Van der Waals",avgDistance=str(round(((interactions[key][1])/interactions[key][3]),3)), avgPotential=str(round(((interactions[key][2])/interactions[key][3]),3))+" kJ/mole", time=str(round(((interactions[key][3]*100)/len(pdbNames)),3))+"%")
			else:
				del interactions[key]

	#removing nodes without adjacent edges
	listToRemove=[]
	for node in graph.nodes():
		if graph.degree[node]==0:
			listToRemove.append(node)
	graph.remove_nodes_from(listToRemove)


	if len(graph.nodes())>1:
			
		#now we will save the graph in a specific format
		if gFormat.lower()=="gml":
			sg.gml(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".gml")
		if gFormat.lower()=="edgelist":
			sg.edgelist(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".edgelist")
		if gFormat.lower()=="gefx":
			sg.gexf(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".gexf")	
		if gFormat.lower()=="graphml":
			sg.graphml(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".graphml")
		if gFormat.lower()=="pajek":
			sg.pajek(graph, folder+"/RIP-MD_Results/Graphs/"+graphName+".pajek")	
	return graph
	
########################################################################
##
## Main function to do the graphs
##
########################################################################
	
def doGraph(pdbDict, outFolder, Format, args, nproc, Time):
	global pdbNames, folder, nodeDict, gFormat, nodeGraph, time
	pdbNames=pdbDict
	folder=outFolder
	gFormat=Format
	time=Time
	############################################
	##
	## Opening node file and generating
	## the node dictionary
	##
	############################################ar
	nodeFile=open(folder+"/RIP-MD_Results/nodes","r")
	nodeDict={}
	for line in nodeFile:
		if line[:-1]!="ID\tname":
			splitted=line[:-1].split("\t")
			nodeDict[splitted[0]]=splitted[1]

	nodeFile.close()
	###############################################
	##
	## And now we will create multigraphs with 
	## the nodes and their properties
	##
	###############################################
	nodeGraph=nx.MultiGraph()
	finalGraph=nx.MultiGraph()
	returnedGraphs=None
	listOfNodesWithAttrs=[]
	if len(pdbNames)==1: # if we have a PDB
		attrFile=open(folder+"/RIP-MD_Results/Attrs/frame_0.attr","r")
		for attr in attrFile:
			splitted=attr[:-1].split("\t") #splitted[0]=nodeID; splitted[1]=sasa; splitted[2]= eight sec structure; splitted[3]= three sec structure
			try:
				listOfNodesWithAttrs.append(nodeDict[splitted[0]])
				nodeData=(nodeDict[splitted[0]]).split(":")
				nodeGraph.add_node(nodeDict[splitted[0]], resname=nodeData[0], chain= nodeData[1], resnum=nodeData[2], eightClasses=splitted[2], threeClasses=splitted[3], SASA=splitted[1])
				finalGraph.add_node(nodeDict[splitted[0]], resname=nodeData[0], chain= nodeData[1], resnum=nodeData[2], eightClasses=splitted[2], threeClasses=splitted[3], SASA=splitted[1])
			except:
				pass	
		for node in nodeDict.items():
			if node[1] not in listOfNodesWithAttrs:
				aux=node[1].split(":")
				nodeGraph.add_node(node[1], resname=aux[0], chain= aux[1], resnum=aux[2])
				finalGraph.add_node(node[1], resname=aux[0], chain= aux[1], resnum=aux[2])
				
		## generating a list with arguments to call the parallel function
		parallelArg=[["(Ca)",args.calpha,"Calpha"], ["(HB)",args.hbond,"HBond"], ["(salt)",args.salt,"Salt"], ["(SS)",args.disulfide,"Disulfide"], ["(cation-pi)",args.cation_pi,"Cation-pi"], ["(pi-pi)",args.pi_pi,"Pi-Pi"], ["(Arg-Arg)",args.arg_arg,"Arg-Arg"], ["(coulomb)",args.coulomb,"Coulomb"], ["(vdW)",args.vdw,"vdW"]]
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		returnedGraphs=list(pool.map(parallelGraph,parallelArg,chunksize=1))
		pool.close()
	
	else:
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		nodeAttrs=list(pool.map(MDProperties,nodeDict.items(),chunksize=1)) #nodeAttrs = [...[node, eightClasses, threeClasses, sasa]...]
		pool.close()

		for attrs in nodeAttrs:
			nodeData=attrs[0].split(":")
			eight1=str(attrs[1])[1:-1]
			aux=eight1.split(",")
			eight=""
			for item in aux:
				if item.split(":")[1]!=" 0.0":
					eight+=(item.split(":"))[0]+":"+(item.split(":"))[1]+"%, "
				

			nodeGraph.add_node(attrs[0],resname=nodeData[0], chain= nodeData[1], resnum=nodeData[2], eightClasses=str(attrs[1])[1:-1], threeClasses=str(attrs[2])[1:-1], avgSASA=str(attrs[3]))
			finalGraph.add_node(attrs[0],resname=nodeData[0], chain= nodeData[1], resnum=nodeData[2], eightClasses=str(attrs[1])[1:-1], threeClasses=str(attrs[2])[1:-1], avgSASA=str(attrs[3]))
		pool=mp.Pool(processes=int(nproc)) #for multiprocessing
		parallelArg=[["(Ca)",args.calpha,"Calpha"], ["(HB)",args.hbond,"HBond"], ["(salt)",args.salt,"Salt"], ["(SS)",args.disulfide,"Disulfide"], ["(cation-pi)",args.cation_pi,"Cation-pi"], ["(pi-pi)",args.pi_pi,"Pi-Pi"], ["(Arg-Arg)",args.arg_arg,"Arg-Arg"], ["(coulomb)",args.coulomb,"Coulomb"], ["(vdW)",args.vdw,"vdW"]]
		returnedGraphs=list(pool.map(MDGraph,parallelArg,chunksize=1))
		pool.close()

	for graph in returnedGraphs:
		
		for  u,v, graphData in graph.edges(data=True):
			edgeList=[]
			attr = dict( (key, value) for key,value in graphData.items())
			edgeList.append((u,v, attr))
			finalGraph.add_edges_from(edgeList)
	
	
	###########################
	## final graph
	###########################
	#now we will save the graph in a specific format
	if gFormat.lower()=="edgelist":
		sg.edgelist(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus.edgelist")
	if gFormat.lower()=="gefx":
		sg.gexf(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus.gexf")	
	if gFormat.lower()=="graphml":
		sg.graphml(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus.graphml")
	if gFormat.lower()=="pajek":
		sg.pajek(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus.pajek")

	sg.saveList(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus_as_list")
	sg.gml(finalGraph, folder+"/RIP-MD_Results/Graphs/consensus.gml") # Always saves GML. Other formats may lose some data. Needed for commDetUI.
	#############################################
	## returning final graph to do covariance
	#############################################
	return finalGraph

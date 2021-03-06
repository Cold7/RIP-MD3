import matplotlib.pyplot as plt
import sys
from glob import glob
import numpy as np

def type_of_interaction(interaction):
	if interaction == "All interactions":
		return "All"
	elif interaction == "C Alpha":
		return "(Ca)"
	elif interaction == "H Bonds":
		return "(HB)"
	elif interaction == "Salt Bridges":
		return "(salt)"
	elif interaction == "Disulfide Bridges":
		return "(SS)"		
	elif interaction == "Cation - pi interaction":
		return "(cation-pi)"		
	elif interaction == "pi - pi interaction":
		return "(pi-pi)"	
	elif interaction == "Arg - Arg interaction":
		return "(Arg-Arg)"
	elif interaction == "VdW contacts":
		return "(vdW)"
	elif interaction == "Coulomb interaction":
		return "(coulomb)"	
	
		
			
def Hist(interaction, folder, original_int):
	edges = glob(folder+"/RIP-MD_Results/Edges/*.edges")
	numbers = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]	

	dict = {}
	for file in edges:
		f = open(file)
		for line in f:
			aux = ""
			try:
				aux = line.split("\t")[2]
			except:
				pass
			if interaction != "All":
				if interaction in aux:
					if aux not in dict:
						dict[aux] = 1
					else:
						dict[aux] += 1
			else:
				if len(aux) > 0:
					if aux not in dict:
						dict[aux] = 1
					else:
						dict[aux] += 1					
		f.close()
	
	
	for item in dict.items():
		percentage = (item[1]*100.0)/len(edges)
		pos = 0
		if 0 <= percentage < 5:
			pos = 0
		if 5 <= percentage < 10:
			pos = 1
		if 10 <= percentage < 15:
			pos = 2
		if 15 <= percentage < 20:
			pos = 3
		if 20 <= percentage < 25:
			pos = 4
		if 25 <= percentage < 30:
			pos = 5
		if 30 <= percentage < 35:
			pos = 6
		if 35 <= percentage < 40:
			pos = 7
		if 40 <= percentage < 45:
			pos = 8
		if 45 <= percentage < 50:
			pos = 9
		if 50 <= percentage < 55:
			pos = 10
		if 55 <= percentage < 60:
			pos = 11
		if 60 <= percentage < 65:
			pos = 12
		if 65 <= percentage < 70:
			pos = 13
		if 70 <= percentage < 75:
			pos = 14
		if 75 <= percentage < 80:
			pos = 15
		if 80 <= percentage < 85:
			pos = 16
		if 85 <= percentage < 90:
			pos = 17
		if 90 <= percentage < 95:
			pos = 18
		if 95 <= percentage <= 100:
			pos = 19
		numbers[pos] += 1

	my_xticks = ['0-5','05-10','10-15','15-20', '20-25', '25-30', '30-35','35-40','40-45', '45-50', '50-55', '55-60', '60-65','65-70','70-75', '75-80', '80-85','85-90','90-95','95-100' ]
	plt.x_range = my_xticks
	plt.bar(my_xticks, numbers)
	plt.title("Number of "+original_int+" present at certain percentaje of frames")
	plt.xlabel("Percentaje of frames")
	plt.ylabel("Frequency")

#	plt.plot(numbers)
#	plt.xlabel("Percentaje of frames")
#	plt.ylabel("Frequency")
#	plt.xticks(my_xticks)
#	plt.title("Number of "+original_int+" present at certain percentaje of frames")
	
	plt.show()


if __name__ == "__main__":
	
	interaction=sys.argv[1]
	original_int = interaction
	interaction = type_of_interaction(interaction)
	folder=sys.argv[2]
	Hist(interaction, folder, original_int)

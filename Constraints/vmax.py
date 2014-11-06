"""
	Given the protein counts, the genes involved in reactions,
	their kcat values, and the Rosetta Stone which converts
	gene names to proteins accession numbers, then we can 
	determine the upper bound of certain flux reactions. 
	The output is a dictionary in the form of {reaction:
	upper bound}.
"""
import cPickle
import numpy
import re
import math

"""
   This pickle contains a list of all the proteins as listed
   by their accession numbers along with the average protein
   counts and standard deviations after the APEX pipeline.
"""

f = open("prot_data_format_apex.pickle", "rb")
load = cPickle.load(f)
f.close()

#Decode the pickle.
average = load[0]
standardDeviation = load[1]
reference  = load[2]

#In the format {YP_number:[[Averages],[Standard Deviations]]}
proteinCounts = {}

#Allocate the pickle into the dictionary.
x = 0 
while x < len(reference):
	thisAverage = numpy.array(average[x]).tolist()
	thisSD = numpy.array(standardDeviation[x]).tolist()
	proteinCounts[reference[x]] = [thisAverage, thisSD]
	x+=1

"""
   File that contains all the reactions that have turnover values.
   Includes reaction abbreviation, reaction name, gene association,
   E.C. number, and turnover number.
"""

f = open("ListOfEnzymes.txt", "rU")
read = f.readlines()
f.close()

#In the format {Reaction Abbrevation: Gene Association}
geneAssociations = {}

#In the format {Reaction Abbrevation: Kcat Value}
kCat = {}

#Sort the file into the dictionaries.
for line in read:
	tabs = line.split('\t')
	reaction = tabs[0]
	genes = tabs[2]
	geneAssociations[reaction] = genes
	kCat[reaction] = tabs[4]

"""
   A semi-manually curated dictionary which has the Blattner
   Lab Number for genes and their corresponding YP Accession
   Number. In the format {b-number:YP_number}
"""

f = open("RosettaStone.pickle", "rb")
rosettaStone = cPickle.load(f)
f.close()

for timepoint in range(0,9):

	#In the format {b-numbers: Average Count}
	geneCounts = {}

	#Based on the rosseta stone, sort the protein counts average values into the dictionary.
	for bGene in rosettaStone:
		proteinID = rosettaStone[bGene]
		try:
			count = proteinCounts[proteinID][0][timepoint]
			if math.isnan(count) != True:
				geneCounts[bGene] = count
		except KeyError:
			print "Some proteins were not found."
		except IndexError:
			print "What the hell kind of timepoint is that?"

	#In the format {Reaction Abbrevation: Enzyme Counts}
	enzymeCounts = {}

	"""
	   Based on the following logic, decide how many molecules of 
	   enzyme existed for each reaction. If there are multiple proteins
	   necessary to catalyze the reaction, then the lowest value from 
	   that set was taken as the count. If there are multiple proteins
	   that indepedently catalyze the reaction, then the sum of that set
	   was taken as the count.
	"""

	for reaction in geneAssociations:
		genes = geneAssociations[reaction]
		if genes.find('or') < 0 and len(genes) > 8:
			genes = '[' + genes + ']'
		genes = genes.replace('(', '[')
		genes = genes.replace(')', ']')
		genes = genes.replace(' and ', ",")
		genes = genes.replace(' or ', ",")
		genes = re.sub(r'(\w+)', r'"\1"' , genes)
		genes = eval(genes)
		if isinstance(genes, str):
			try:
				count = geneCounts[genes]
				enzymeCounts[reaction] = count
			except:
				#print "Could not find a value for %s gene." % (genes)
				pass
		if isinstance(genes, list):
			count = []
			for ors in genes:
				if isinstance(ors, list):
					andList = []
					for ands in ors:
						try:
							andCount = geneCounts[ands]
							andList.append(andCount)
						except:
							#print "Could not find a value for %s gene." % (ands)
							pass
					try:
						thisCount = min(andList)
						count.append(thisCount)
					except:
						pass
				else:
					try:
						orCount = geneCounts[ors]
						count.append(orCount)
					except:
						#print "Could not find a value for %s gene." % (ors)
						pass

	#In the format {Reaction : Upper Bound}
	vmax = {}

	for reaction in enzymeCounts:
		count = enzymeCounts[reaction]
		thisKCat = float(kCat[reaction])

		avagadro = 6.02214129*(10**23)
		weight = 2.58 * (10**-13)
		cellsPerGram = 1/weight

		bound = count * (1/avagadro) * (10**3) * thisKCat * cellsPerGram * (60**2)

		vmax[reaction] = bound

	fileName = 'Vmax/VMaxTimepoint' + str(timepoint+1) + '.pickle'

	output = open(fileName, 'wb')
	cPickle.dump(vmax, output)
	output.close()




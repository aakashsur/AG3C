
"""
	Given the protein counts, the genes involved in reactions,
	their kcat values, and the Rosetta Stone which converts
	gene names to proteins accession numbers, then we can 
	determine the upper bound of certain flux reactions. 
	The output is a dictionary in the form of {reaction:
	upper bound}.
"""
import pickle
import numpy
import re
import itertools
import os

# Returns a dictionary of kCat values and gene associations for each reaction.
def getInfo ():

	"""
	   File that contains all the reactions that have turnover values.
	   Includes reaction abbreviation, reaction name, gene association,
	   E.C. number, and turnover number.
	"""

	inFile = open("ListOfEnzymes.txt", "rU")
	read = inFile.readlines()
	inFile.close()

	# In the format {Reaction Abbrevation: kCat Value}
	kCat = {}
	# In the format {Reaction Abbrevation: Gene Association}
	geneAssociations = {}

	# Read the tab delimited file and sort the columns.
	i = 0
	for line in read:
		tabs = line.split('\t')
		reaction = tabs[0]
		kCat[reaction] = tabs[4]
		geneAssociations[reaction] = tabs[2]
		i += 1

	print ('There are %i reactions being parsed.') % (i)

	# Parse the gene logic structure. 
	for reaction in geneAssociations:
		genes = geneAssociations[reaction]
		# If it is AND logic with multiple genes. 
		if genes.find('or') < 0 and len(genes) > 8:
			genes = '[' + genes + ']'
		# Make it into list format.
		genes = genes.replace('(', '[').replace(')', ']')
		genes = genes.replace(' and ', ",").replace(' or ', ",")
		# Add quote around all the gene names.
		genes = re.sub(r'(\w+)', r'"\1"' , genes)
		# Evaluate it as a list object.
		genes = eval(genes)
		geneAssociations[reaction] = genes

	return geneAssociations, kCat

"""
	Open the protein list for each confidence level at each timepoint. Since the pickle
	files are in the form of an ascending sorted list of tuples, we need to reverse the 
	mean counts file and then convert both the mean counts and rsd into dictionaries so 
	they can be more easily utilized. We also take the top n proteins as specified by the 
	size variable. 
"""

def getProteinCounts (percent, time, size):

	# Open ranked average protien counts pickle.  
	fileName = ('ProteinList/{0}Confidence/Mean/Time{1}.pickle').format(percent, time+1)
	inFile = open(fileName, 'rb')
	proteinCounts = pickle.load(inFile)[::-1]
	inFile.close()

	print ('Loading data for %i proteins.') % (len (proteinCounts))

	# Open relative standard deviation rankings pickle. 
	fileName = ('ProteinList/{0}Confidence/RSD/Time{1}.pickle').format(percent, time+1)
	inFile = open(fileName, 'rb')
	rsd = pickle.load(inFile)
	inFile.close()

	# Since we want to test the top n values, we also want to see what happens when we include all proteins. 
	if not size:
		size = len(proteinCounts)

	# The top n proteins and their counts.
	# In the format {YP_number: Count}. 
	meanCount = {}
	proteinCountsComplete = {}

	i = 0
	for pair in proteinCounts:
		if i < size:
			meanCount[pair[0]] = pair[1]
		proteinCountsComplete[pair[0]] = pair[1]
		i += 1

	# The top n rsd values to create a dictionary. 
	# In the format {YP_number: Count}.
	rsdCount = {}
	for pair in rsd[:size]:
		rsdCount[pair[0]] = proteinCountsComplete[pair[0]]

	return proteinCountsComplete, meanCount, rsdCount

def getGeneCounts (rosettaStone, proteinCounts, verbose = False):

	# In the format {b-numbers: Average Count}
	geneCounts = {}

	# Based on the rosetta stone, sort the protein counts average values into the dictionary.
	for gene in rosettaStone:
		proteinID = rosettaStone[gene]
		if proteinID in proteinCounts:
			geneCounts[gene] = proteinCounts[proteinID]
		elif verbose:
			print ("%s not found in protein counts set.") % (gene)

	print ('There are %i proteins which are not in Rosetta Stone but have counts.') % (len (set (proteinCounts) - set (rosettaStone.values())))
	print ('There are %i genes for which there is currently count data.\n') % (len (geneCounts))

	return geneCounts

"""
   Based on the following logic, decide how many molecules of 
   enzyme existed for each reaction. If there are multiple proteins
   necessary to catalyze the reaction, then the lowest value from 
   that set was taken as the count. If there are multiple proteins
   that indepedently catalyze the reaction, then the sum of that set
   was taken as the count.
"""

def getEnzymeCounts (rosettaStone, proteinCountsComplete, proteinCounts, geneAssociations, geneCounts, info, verbose = False):

	#In the format {Reaction Abbrevation: Enzyme Counts}
	enzymeCounts = {}

	for reaction in geneAssociations:
		genes = geneAssociations[reaction]
		# If there is only one gene for the reaction.
		if isinstance(genes, str):
			if genes in geneCounts:
				enzymeCounts[reaction] = geneCounts[genes]
			elif verbose and genes not in rosettaStone:
				print ("Gene %s is not in Rosetta Stone.") % (genes)
			elif verbose and rosettaStone[genes] not in proteinCounts:
				print ('Protein %s is not in this set.') % (rosettaStone[genes])

		# If there is multiple genes involved, parse the gene logic. 
		elif isinstance(genes, list):
			# List of all possible enzyme units that can catalyze the reaction.
			# May be multiple genes in the unit. 
			count = []
			# First level, only OR logic, i.e. multiple alternative enzymes. 
			for ors in genes:
				# If there is a second level of logic. 
				if isinstance(ors, list):
					andList = []
					# Second level, only AND logic, i.e. multiple subunits of the enzyme. 
					for ands in ors:
						if ands in geneCounts:
							andList.append(geneCounts[ands])
						elif verbose and ands not in rosettaStone:
							print ("Gene %s is not in Rosetta Stone.") % (ands)
						elif verbose and rosettaStone[ands] not in proteinCounts:
							print ('Protein %s is not in this confidence level.') % (rosettaStone[ands])

					if len (andList) > 0: 
						count.append(min(andList))

				elif info.lower() == 'mean':
					if (ors in rosettaStone) and (rosettaStone[ors] in proteinCountsComplete):
						count.append(proteinCountsComplete[rosettaStone[ors]])
					elif verbose and ors not in rosettaStone:
						print ("Gene %s is not in Rosetta Stone.") % (ors)
					elif verbose and rosettaStone[ors] not in proteinCountsComplete:
						print ('Protein %s is not in this confidence level.') % (rosettaStone[ors])

				else:
					if ors in geneCounts:
						count.append(geneCounts[ors])
					elif verbose and ors not in rosettaStone:
						print ("Gene %s is not in Rosetta Stone.") % (ors)
					elif verbose and rosettaStone[ors] not in proteinCounts:
						print ('Protein %s is not in this set.') % (rosettaStone[ors])

			if sum (count) > 0:
				enzymeCounts[reaction] = sum(count)

		else:
			print ('Failure in gene association logic.')

	print ('There are %i reactions which have enzyme counts.') % (len (enzymeCounts))
	return enzymeCounts

def getVMax (enzymeCounts, kCat, percent, time, size, info, write = False):
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

	if write:
		if not size:
			size = 'Full'
		folderName = ('Vmax/{0}Confidence/{1}/Size{2}').format(percent, info, size)
		if not os.path.isdir(folderName):
			os.makedirs(folderName)
		fileName = folderName + ('/Time{}.pickle').format(time + 1)
		# Don't overwrite files. 
		outFile = open(fileName, 'wb')
		pickle.dump(vmax, outFile)
		outFile.close()

	print ('There are %i reactions which have upper bounds.') % (len (enzymeCounts))
	return vmax

def main ():

	"""
	   A semi-manually curated dictionary which has the Blattner
	   Lab Number for genes and their corresponding YP Accession
	   Number. In the format {b-number:YP_number}.   
	   Use the pickle with a list of all the proteins as listed
	   by their accession numbers along with the average protein
	   counts and standard deviations after the APEX pipeline.
	"""

	inFile = open("RosettaStone.pickle", "rb")
	rosettaStone = pickle.load(inFile)
	inFile.close()

	geneAssociations, kCat = getInfo()

	# Iterate through each confidence level.
	confidence = [90, 95, 99]
	# Iterate through different list sizes.
	limit = [50, 100, 200, False]
	# Iterate through each timepoint.
	timepoint = range (9)

	for percent, size, time in itertools.product(confidence, limit, timepoint):
		print (percent, time + 1, size)

		proteinCountsComplete, meanCount, rsdCount = getProteinCounts(percent, time, size)

		geneMean = getGeneCounts(rosettaStone, meanCount)
		geneRSD = getGeneCounts(rosettaStone, rsdCount)

		enzymeMean = getEnzymeCounts(rosettaStone, proteinCountsComplete, meanCount, geneAssociations, geneMean, 'mean')
		enzymeRSD = getEnzymeCounts(rosettaStone, proteinCountsComplete, rsdCount, geneAssociations, geneRSD, 'rsd')

		vmaxMean = getVMax(enzymeMean, kCat, percent, time, size, 'Mean', True)
		vmaxRSD = getVMax(enzymeRSD, kCat, percent, time, size, 'RSD', True)

main ()



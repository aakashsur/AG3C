"""
	We need to assign protein counts data, but all we have are 
	the b-number associated with each reaction. To connect the
	dots, we ask Entrez what gene symbols are associated with 
	the b-number and what the protein names of those gene symbols 
	are. The next script, rosettaStone.py will ask Entrez questions 
	about the protein accession number and use the information 
	generated from this program to assign proteins to b-numbers.
"""
from Bio import Entrez
import time
import re
import cPickle

"""
   File that contains all the reactions that have turnover values.
   Includes reaction abbreviation, reaction name, gene association,
   E.C. number, and turnover number.
"""

f = open("ListOfEnzymes.txt", "rU")
read = f.readlines()
f.close()

#Get all Blattner Lab Number Id's for the reactions.

regex = r'b\d\d\d\d' 
listOfGenes={} 

for line in read:
	tabs = line.split('\t')
	genes = tabs[2]
	genes = re.findall(regex, genes)
	for gene in genes:
		listOfGenes[gene] = []

Entrez.email = "aakash.sur@utexas.edu" #For Enztrez to contact you if your querys are causing problems. 

"""
   This loop will go through each Blattner Number Gene and search 
   the Gene database on NCBI's website. The results are then 
   passed to a second loop which will get more detailed information 
   about each search result. The idea is to gather information about 
   each b-number that will help us identify how the gene relates to
   the protein accession number. The path is b-number -> gene symbol
   -> protein accession. This script will gather the gene name and
   its synonyms, and the corresponding protein names and E.C. numbers.
   rosettaStone.py will use this information to connect the dots.
""" 

for gene in listOfGenes: #For each b-number gene:

	print gene

	thisGene = gene+'[Gene]'
	handle = Entrez.esearch(db="gene", term= thisGene)
	record = Entrez.read(handle)
	idList = record["IdList"]

	synonymList = []
	proteinNameList = []
	proteinECList = []

	for idNumber in idList:

		time.sleep(1) #NCBI wants to limit access to one query per second. (Biopython will enforce this.)
		handle = Entrez.efetch(db="gene", id= idNumber, retmode='xml')
		record = Entrez.read(handle)

		#Since we don't know whether each record will have all the parts, I try to catch all the errors.
		try:
			geneName = record[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
		except:
			geneName = "No gene name available."
		try:
			synonym = record[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_syn"]
		except:
			synonym = "No synonyms available."
		try:
			proteinName = record[0]["Entrezgene_prot"]['Prot-ref']['Prot-ref_name'][0]
		except:
			proteinName = "No protein name available."
		try:
			proteinEC = record[0]["Entrezgene_prot"]['Prot-ref']['Prot-ref_ec'][0]
		except:
			proteinEC = "No E.C. available."

		geneNameList.append(geneName)

		#We are only interested in the three to four letter gene symbols, not gene IDs.
		for name in synonym:
			if "JW" in name:
				pass
			elif "ECK" in name:
				pass
			elif name in gene:
				pass
			else:
				geneNameList.append(name)

		proteinNameList.append(proteinName)
		proteinECList.append(proteinEC)

		#We only want unique names in our dictionary, and if there is information available, we will disregard
		#the fact that some of the search queries for a b-number are missing information. 
		if geneNameList.count("No gene name available.") == len(idList):
			geneNameList = "No gene name available."
		else:
			geneNameList = [item for item in geneNameList if item != "No gene name available."]
			geneNameList = list(set(geneNameList))

		if proteinNameList.count("No protein name available.") == len(idList):
			proteinNameList = "No protein name available."
		else:
			proteinNameList = [item for item in proteinNameList if item != "No protein name available."]
			proteinNameList = list(set(proteinNameList))

		if proteinECList.count("No E.C. available.") == len(idList):
			proteinECList = "No E.C. available."
		else:
			proteinECList = [item for item in proteinECList if item != "No E.C. available."]
			proteinECList = list(set(proteinECList))

	listOfGenes[gene] = [geneNameList, proteinNameList, proteinECList]

#In the format {b-number: [[gene name and synonyms], [protein names], [E.C. numbers]]}
print listOfGenes

#Yum, pickles. Since this searching Entrez takes time, let's pickle the results!
output = open('geneInfoList.pickle', 'wb')
cPickle.dump(listOfGenes, output)
output.close()





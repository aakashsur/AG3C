import cPickle
import numpy
from Bio import Entrez
import time

"""
   This is a list of all the RefSeq Accession Numbers and 
   their corresponding gene symbols and MRN ID numbers.
"""

f = open('REL606_gene_ids.txt', 'rU')
read = f.readlines()
f.close()

geneConversionList = {} 
#In the format {Gene Symbol: YP_accession Number}

for line in read:
	tabs = line.split('\t')
	geneConversionList[tabs[1]]= tabs[0]

"""
   File that contains all the reactions that have turnover values.
   Includes reaction abbreviation, reaction name, gene association,
   E.C. number, and turnover number.
"""

f = open("ListOfEnzymes.txt", "rU")
reactions = f.readlines()
f.close()

"""
   This is a dictionary of all the Blattner Number Genes
   involved in reactions with turnovers. It has information
   about gene names, gene synonyms, protein names, and E.C.
   numbers.

   It contains {b-number:[Gene Names],[Protein Names],[E.C. Numbers]}
"""

f = open("geneInfoList.pickle", "rb")
geneInfoList = cPickle.load(f)
f.close()

Entrez.email = "aakash.sur@utexas.edu"

rosettaStone = {} 
#In the format {b-number:YP_accession Number}

matchDict = {} 
#In the format {b-number:[Gene Symbols]}

"""
   The primary matching portion of the script. The loop takes all
   the Blattner Number genes from the geneInfoList.pickle and will
   go through each of its gene names and synonyms, and will try and
   match it to a protein via the REL606_gene_ids.txt. However, trouble
   is always near by: sometimes there are b-numbers with no matches,
   and sometimes there are b-numbers with more than one match. At the
   moment, I don't do anything for non-matching b-numbers, but I do 
   catch all the genes that map to multiple proteins. We will sort this
   mess out below.
"""

for bGene in geneInfoList:
	geneTestList = geneInfoList[bGene][0]
	matchDict[bGene] = []
	matchNumber = 0

	for geneSymbol in geneTestList:

		try:
			rosettaStone[bGene] = geneConversionList[geneSymbol]
			matchDict[bGene].append(geneSymbol)
			matchNumber += 1
		except:
			pass

	if matchNumber == 0:
		print bGene + " has no matches."

print "Total number of gene and proteins matched is " + str(len(rosettaStone))

geneChoices = {} 
#In the format {b-number: [[[YP_number, Protein Name, E.C. Number] ,[2nd Protein]], [Reactions]]}

"""
   So the correct protein for a gene is often apparent once we know
   the name of that protein. The idea here is to search Entrez for
   the YP_accession number, which should only return one result,
   then fetch the details of that result and pull out the protein 
   name and E.C. number. Then, we can pull out the reactions that
   are related to the gene in question and manually go through and
   assign which protein ought to map to which gene.
"""

for bGene in rosettaStone:
	geneMatchList = matchDict[bGene]

	if len(geneMatchList) > 1:
		geneChoices[bGene] = [[]]

		for geneSymbol in geneMatchList:
			proteinID = geneConversionList[geneSymbol]

			handle = Entrez.esearch(db="protein", term= proteinID)
			record = Entrez.read(handle)
			idList = record["IdList"]
			idNumber = idList[0]

			time.sleep(1) #NCBI wants to limit access to one query per second. (Biopython will enforce this.)
			handle = Entrez.efetch(db="protein", id= idNumber, retmode='xml')
			record = Entrez.read(handle)

			try:
				proteinName = record[0]['GBSeq_feature-table'][1]['GBFeature_quals'][0]['GBQualifier_value']
			except:
				proteinName = "No protein name."
			try:
				ecNumber = record[0]['GBSeq_feature-table'][1]['GBFeature_quals'][1]['GBQualifier_value']
			except:
				ecNumber = "No E.C. Number."

			proteinList = [proteinID, proteinName, ecNumber]
			geneChoices[bGene][0].append(proteinList)

		reactionList = [line.split('\t')[1] + ' E.C. ' + line.split('\t')[3] for line in reactions if bGene in line]
		geneChoices[bGene].append(reactionList)

		print "Checked: " + bGene

print geneChoices

"""
   Nobody wants to do things by hand. Lets build an user interface
   that will help us do this. Here, for each b-number gene that
   mapped to multiple proteins, the reaction is first displayed,
   then the protein choices. The user is then prompted to choose
   which protein they would like to be assigned to that gene.
"""

for bGene in geneChoices:
	print "The reactions are: "
	for line in geneChoices[bGene][1]:
		print line
	count = 1

	for protein in geneChoices[bGene][0]:
		thisOption = "Option " + str(count) + ": " + str(protein)
		print thisOption
		count += 1

	print 'Which option would you like?'
	choice = input("Enter Number ")
	choice -= 1
	ypNumber = geneChoices[bGene][0][choice][0]
	rosettaStone[bGene] = ypNumber

	print "Now " + bGene + " has a value of " + ypNumber + "\n"

print rosettaStone

"""
   Since there are time intensive Entrez searches and a manual curation
   step, nobody wants to wait around. So let's have a pickle party.
   Pickle everything!
"""

output = open('RosettaStone.pickle', 'wb')
cPickle.dump(rosettaStone, output)
output.close()

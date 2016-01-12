
from Bio import Entrez
import re
import pickle
import time
import os

def getGeneInfo (geneID):

	# Search the protein database for gene ID.
	handle = Entrez.esearch(db='protein', term = geneID)
	record = Entrez.read(handle)

	# If there are multiple results throw an error, otherwise set the ID. 
	if len(record['IdList']) > 1:
		print (('More than one record for: {}').format(geneID))
		return
	else:
		ID = (record["IdList"][0])

	# Now find the actual entry in the protein database with all the details. 
	handle = Entrez.efetch(db='protein', id = ID, retmode = 'xml')
	handle = Entrez.read (handle)

	# Get the portion about the Blattner number. 
	regex = r'b\d\d\d\d' 
	genes = re.findall(regex, str(handle))
	if len (genes) == 1:
		return genes[0]
	elif len (genes) > 1:
		print ('More than one b-number for: {}').format(geneID)
	else:
		print ('Does not have a b-number: {}').format(geneID)

def main ():

	Entrez.email = 'aakash.sur@utexas.edu'
	inFile = open ('REL606_gene_ids.txt', 'rU')
	inFile.readline()

	outFile = open('RosettaStone.txt', 'a+')
	outFile.seek(0)
	read = outFile.readlines()

	if len (read) == 0:
		outFile.write ('protein ID\tmRNA ID\tBlattner Number\n')
		lastProtein = None
		valid = True
	else:
		lastProtein = read[-1].split('\t')[0]
		valid = False

	for line in inFile:
		line = line.strip().split('\t')
		protein = line[0]
		if lastProtein == protein:
			valid = True
		if not valid:
			continue 
		gene = line [2]
		answer = getGeneInfo (gene)
		string = ('{}\t{}\t{}\n').format (protein, gene, answer)
		outFile.write (string)
		print (string.strip())
		time.sleep (0.35)

	inFile.close()
	outFile.close()

main ()


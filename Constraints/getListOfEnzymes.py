"""
	A list of the reactions from the iJO1366 Model created
	by the Paulson lab. The list contains reaction abbrevations, 
	reaction names, gene associations, E.C. numbers, turnover
	rates, substrate information, species, commentary, and 
	reference. This was downloaded from a paper by Labhsetwar, 
	et al. titled Heterogeneity in protein expression induces 
	metabolic variability in a modeled Escherichia coli population - 
	which explores using protein levels to constrain an FBA.
"""

f = open("ReactionsWithTurnover.txt", "rU")
read = f.readlines()
f.close()

turnoverDict = {}

for line in read:
	tabs = line.split('\t')
	if tabs[4]:
		turnoverDict[tabs[0]] = [tabs[1], tabs[2], tabs[3], tabs[4], tabs[5]]


"""
   Create a new tab-delimted file containing only reaction
   information for those reactions that have an emperical
   turnover rate. 
"""

f = open("ListOfEnzymes.txt", "w")

for item in turnoverDict:
	line = item
	for info in turnoverDict[item]:
		line += "\t" + info 
	line += "\n"
	f.write(line)

f.close()


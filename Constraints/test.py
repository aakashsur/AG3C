
inFile = open('RosettaStone.txt')

inFile.readline()
out = {}

for line in inFile:
	line = line.strip().split('\t')
	gene = line[2].split(',')
	if gene[0] == 'None':
		continue
	protein = line[0]
	for item in gene:
		out[item] = protein

print out
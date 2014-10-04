import cPickle
import numpy

f = open("VMax.pickle", "rb")
vMax = cPickle.load(f)
f.close()

print vMax['PRFGS']

f = open("RosettaStone.pickle", "rb")
rosettaStone = cPickle.load(f)
f.close()

print rosettaStone['b2557']

f = open("prot_data_format_apex.pickle", "rb")
load = cPickle.load(f)
f.close()

average = load[0]
standardDeviation = load[1]
reference  = load[2]

proteinCounts = {}
#In the format {YP_number:[[Averages],[Standard Deviations]]}

x = 0 

while x < len(reference):
	thisAverage = numpy.array(average[x]).tolist()
	thisSD = numpy.array(standardDeviation[x]).tolist()
	proteinCounts[reference[x]] = [thisAverage, thisSD]
	x+=1


print proteinCounts['YP_003045639.1']
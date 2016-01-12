import cPickle
import numpy as np
import matplotlib.pyplot as plt

# Iterate through each confidence level.
for percent in [90,95,99]:
	print (percent)

	# Iterate through each timepoint. 
	for time in xrange (9):

		fileName = '../Constraints/ProteinList/'+ str(percent) + 'Confidence/' + 'Mean/Time' + str(time+1) + '.pickle'
		f = open(fileName, 'rb')
		meanCounts = cPickle.load(f)
		f.close()

		fileName = '../Constraints/ProteinList/'+ str(percent) + 'Confidence/' + 'RSD/Time' + str(time+1) + '.pickle'
		f = open(fileName, 'rb')
		rsd = cPickle.load(f)
		f.close()

		counts = {}
		for protein in meanCounts:
			counts[protein[0]] = protein[1]

		x = []
		y = []

		for protein in rsd:
			x.append(protein[1])
			y.append(counts[protein[0]])

		print ('There are %i proteins in the %i timepoint.') % (len(meanCounts), time + 1)

		plt.scatter(x, y)
		plt.show()




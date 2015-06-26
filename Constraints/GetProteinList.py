
import cPickle
import numpy as np
import operator
import os

def writeOufFile (percent, time, data, info):

	# Save it as a pickle file. 
	folderName = ('ProteinList/{}Confidence/{}').format(percent, info)
	if not os.path.isdir(folderName):
		os.makedirs(folderName)
	fileName = folderName + ('/Time{}.pickle').format(time + 1)
	if not os.path.isfile(fileName):
		outFile = open(fileName, 'wb')
		cPickle.dump(data, outFile)
		outFile.close()

def main ():

	# Pickle containing all the protein counts information. 
	# In the format [[Average Counts], [Standard Deviations], [YP_number]]
	inFile = open("prot_data_format_apex.pickle", "rb")
	data = cPickle.load(inFile)
	inFile.close()

	# 90% Confidence = 1.2815, 95% Confidence = 1.645, 99% Confidence = 2.327
	level = {90: 1.2815, 95: 1.645, 99 : 2.327}

	# Iterate through each confidence level.
	for percent in level:
		print (percent)
		zScore = level[percent]
		# Iterate through each timepoint. 
		for time in xrange (9):
			print (time)
			# Dictionary in the format {YP_number:Average} for the timepoint. 
			mean = {}
			# Dictionary in the format {YP_number:Relative Standard Deviation} for the timepoint. 
			rsd = {}
			# Iterate through each protein.
			for x in xrange (len (data[0])):
				average = data[0][x]
				sd = data[1][x]
				name = data[2][x]

				# Toss out the protein if there are any unspecified values or if it fails the one-tailed test. 
				if np.isnan(average).any() or ((average[time] - zScore*sd[time]) <= 0):
					continue

				mean[name] = average[time]
				# Calculate the relative standard deviation.
				rsd[name] = (sd[time]*100/average[time])

			# Sort the dictionaries into a list of tuples going in ascending order sorted by the values.
			# In the format [('YP_number', count)]
			mean = sorted (mean.items(), key = operator.itemgetter (1))
			rsd = sorted (rsd.items(), key = operator.itemgetter (1))

			# Write results into a pickle
			writeOufFile (percent, time, mean, 'Mean')
			writeOufFile (percent, time, mean, 'RSD')

main()



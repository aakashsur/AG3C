from matplotlib import pylab as plt
import numpy as np
import pandas as pd
import itertools

confidence = [90,95,99]
info = ['Mean', 'RSD']
size = [50,100,200,'Full']

def main():
	for percent, kind, num in itertools.product(confidence, info, size):
		# print (percent, kind, num)
		fileName = ('../Constraints/Results/{}Confidence/{}/Size{}.txt').format(percent, kind, num)
		title = ('{} Confidence {} Size {}').format(percent, kind, num)
		outFile = ('TimeSeries/{}Confidence{}Size{}').format(percent, kind, num)

		df = pd.read_table(fileName)
		growth = df['Growth Rate']
		data = df.drop(['Growth Rate'], axis = 1)

		data.plot()
		fig = plt.figure(1)
		ax = fig.add_subplot(111)
		plt.xlabel('Timepoint')
		plt.ylabel('mmol/g*hr')
		plt.title(title)
		handles, labels = ax.get_legend_handles_labels()
		lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.2,1))
		fig.savefig(outFile, bbox_extra_artists=(lgd,), bbox_inches='tight')
		plt.close()

main()
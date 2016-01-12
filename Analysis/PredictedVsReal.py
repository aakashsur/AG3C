from matplotlib import pylab as plt
import numpy as np
import pandas as pd
import itertools

database = []

confidence = [90,95,99]
info = ['Mean', 'RSD']
size = [50,100,200,'Full']

def main():
	x = pd.read_table('RealFlux.txt')

	for percent, kind, num in itertools.product(confidence, info, size):
		# print (percent, kind, num)
		fileName = ('../Constraints/Results/{}Confidence/{}/Size{}.txt').format(percent, kind, num)


		
		df = pd.read_table(fileName)
		y = df.drop(['Growth Rate'], axis = 1)

		n = 0 
		colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
		for column in x:
			xAxis = x[column].values
			yAxis = y[column].values
			plt.scatter(xAxis,yAxis, c = colors[n], label = column)
			m,b = np.polyfit(xAxis, yAxis, 1) 
			plt.plot(xAxis, m*xAxis + b, c = colors[n])
			n += 1

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
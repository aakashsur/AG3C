from matplotlib import pylab as plt
import numpy as np

f = open('FluxData.txt', 'rU')
data = {}
x = 1
for line in f:
	if x == 9:
		x = 1
	if x == 1:
		header = line.strip()
		data[header] = {}
		x+=1
		continue
	line = line.strip()
	line = line.split('\t')
	data[header][line[0]] = map(float, line[1:])
	x+=1

x = [1,2,3,4,5,6,7,8,9]

# correlations = {}

# for reaction in data['Real']:
# 	correlations[reaction] = {}
# 	for condition in data:
# 		correlations[reaction][condition] = data[condition][reaction]
# colors = ['b', 'g', 'r', 'c', 'm', 'y', 'b']
# f = open('correlations.txt' , 'w')
# for reaction in correlations:
# 	f.write(reaction + '\n')
# 	i = 0
# 	for condition in correlations[reaction]:
# 		x = correlations[reaction]['Real']
# 		y = correlations[reaction][condition]

# 		plt.scatter(x, y, c = colors[i], label = condition)
# 		m,b = plt.polyfit(x, y, 1) 
# 		yp = plt.polyval([m,b],x)
# 		plt.plot(x, yp, c = colors[i])

# 		f.write(condition + '\t' +str(m) + '\n')

# 		i+=1
# 	fig = plt.figure(1)
# 	ax = fig.add_subplot(111)
# 	plt.xlabel('Real Flux')
# 	plt.ylabel('Predicted Flux')
# 	plt.title('Comparison of ' + reaction)
# 	handles, labels = ax.get_legend_handles_labels()
# 	lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.4,1))
# 	fig.savefig(reaction, bbox_extra_artists=(lgd,), bbox_inches='tight')
# 	plt.close()


# for condition in data:
# 	for reaction in data[condition]:
# 		plt.plot(x,data[condition][reaction], label = reaction)
# 	fig = plt.figure(1)
# 	ax = fig.add_subplot(111)
# 	plt.xlabel('Timepoint')
# 	plt.ylabel('mmol/g*hr')
# 	plt.title('Flux ' + condition)
# 	handles, labels = ax.get_legend_handles_labels()
# 	lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.2,1))
# 	fig.savefig(condition, bbox_extra_artists=(lgd,), bbox_inches='tight')
# 	plt.close()






import pickle
import numpy as np

# Pickle containing all the protein counts information. 
# In the format [[Average Counts], [Standard Deviations], [YP_number]]
f = open("prot_data_format_apex.pickle", "rb")
data = pickle.load(f)
f.close()

count = 0
for item in data[0]:
	if np.isnan(item[:7]).any():
		count += 1

print (len (data[0]))
print (count)


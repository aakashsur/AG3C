import cPickle
from cobra.io import *
from cobra.flux_analysis import *

f = open("VMax.pickle", "rb")
vMax = cPickle.load(f)
f.close()

model = read_sbml_model("REL.xml")
model.optimize()
print model.solution.f

for reaction in vMax:
	upperBound = vMax[reaction] * 100
	try:
		rxn = model.reactions.get_by_id(reaction)
		rxn.upper_bound = upperBound
	except KeyError: 
		print reaction

model.optimize()
print model.solution.f

# for rxn in model.reactions:
# 	if "EX_" in str(rxn):
# 		print rxn

# # reaction = 'EX_glc_LPAREN_e_RPAREN_'
# # rxn = model.reactions.get_by_id(reaction)
# # print rxn
# # rxn.upper_bound = 99999
# # rxn.lower_bound = -99999

# model.optimize()
# print model.solution.f

    
	
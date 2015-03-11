"""
	The last script in the sequence. Adds the bounds specified
	in the pickle files in the folder /Vmax. Output two text
	files which have the flux ratios for the specified reactions
	over all the time points, one for the unconstrained fluxes 
	one for the constrained fluxes. 
"""
import cPickle
from cobra.io import *
from cobra.flux_analysis import *
from escher import Builder

f = open('RegularFluxes.txt', 'w')

model = read_sbml_model("REL.xml")
model.optimize()


# We lower the default bounds to -1,000 and 1,000.
newBounds = {999999.000000: 1000.0, -999999.0 : -1000.0}

for rxn in model.reactions:
	if "EX_" in str(rxn):
		thisReaction = str(rxn)
		try:
			rxn.upper_bound = newBounds[rxn.upper_bound]
			rxn.lower_bound = newBounds[rxn.lower_bound]
		except KeyError:
			pass
model.optimize()
print model.solution.f

# Equations for flux ratios. 
solution = model.solution.x_dict
ppck = solution['PPCK']
eno = solution['ENO']
PEPfromOAA = ppck/(ppck+eno)
thisLine = "PEP from OAA\t" + str(PEPfromOAA) + "\n"
f.write(thisLine)

tkt2 = solution['TKT2']
tala = solution['TALA']
E4PThroughTK = tkt2/(tkt2+tala)
thisLine = "E4P through TK\t" + str(E4PThroughTK)+ "\n"
f.write(thisLine)

ppc = solution['PPC']
mdh = solution['MDH']
mdh2 = solution['MDH2']
mdh3 = solution['MDH3']
OAAfromPEP = ppc/(ppc+mdh+mdh2+mdh3)
thisLine = "OAA from PEP\t" + str(OAAfromPEP)+ "\n"
f.write(thisLine)

tkt1 = solution['TKT1']
eda = solution['EDA']
pfk = solution['PFK']
PEPthroughTK = tkt1/(eda + pfk + tala + tkt1)
thisLine = "PEP through TK\t" + str(PEPthroughTK)+ "\n"
f.write(thisLine)

ghmt2r = solution['GHMT2r']
serat = solution['SERAT']
lserdhr = solution['LSERDHr']
psp_l = solution['PSP_L']
sert4pp = solution['SERt4pp']
sert2rpp = solution['SERt2rpp']
SERfromGLY = ghmt2r/(ghmt2r+serat+lserdhr+psp_l+sert4pp+sert2rpp)
thisLine = "SER from GLY\t" + str(SERfromGLY)+ "\n"
f.write(thisLine)

glyat = solution['GLYAT']
pragsr = solution['PRAGSr']
amptasecg = solution['AMPTASECG']
amptasepg = solution['AMPTASEPG']
sarcox = solution['SARCOX']
thra2i = solution['THRA2i']
thrai = solution['THRAi']
GLYfromSER = ghmt2r/(ghmt2r+glyat+pragsr+amptasecg+amptasepg+sarcox+thra2i+thrai)
thisLine = "GLY from SER\t" + str(GLYfromSER)+ "\n"
f.write(thisLine)

me1 = solution['ME1']
me2 = solution['ME2']
pyk = solution['PYK']
PYRfromMAL = (me1+me2)/(me1+me2+pyk+eda)
thisLine = "PYR from MAL\t" + str(PYRfromMAL)+ "\n"
f.write(thisLine)
f.close()

lines = ["PEP from OAA\t","E4P through TK\t","OAA from PEP\t","PEP through TK\t","SER from GLY\t","GLY from SER\t", "PYR from MAL\t"]
outLine = []

# Calculating fluxes for the constrained models.
for timepoint in range(0,9):

	print timepoint

	# Tells script where to pull pickle.
	fileName = 'Vmax/VMaxTimepoint' + str(timepoint+1) + '.pickle'

	f = open(fileName, "rb")
	vMax = cPickle.load(f)
	f.close()

	model = read_sbml_model("REL.xml")
	model.optimize()

	# We lower the default bounds to -1,000 and 1,000.
	newBounds = {999999.000000: 1000.0, -999999.0 : -1000.0}

	for rxn in model.reactions:
		if "EX_" in str(rxn):
			thisReaction = str(rxn)
			try:
				rxn.upper_bound = newBounds[rxn.upper_bound]
				rxn.lower_bound = newBounds[rxn.lower_bound]
			except KeyError:
				pass

	# for reaction in vMax:
	# 	upperBound = vMax[reaction]
	# 	try:
	# 		rxn = model.reactions.get_by_id(reaction)
	# 		oldBound = rxn.upper_bound
	# 		rxn.upper_bound = upperBound
	# 		model.optimize()
	# 		print ("Changing %s results in %f")%(rxn, model.solution.f)
	# 		rxn.upper_bound = oldBound
	# 	except KeyError: 
	# 		pass
	# 		#print reaction

	print "Changing bounds."

	for reaction in vMax:
		if reaction == "GMPS2" or reaction == 'PRFGS' or reaction == 'ADSL2r':
			# pass
		upperBound = vMax[reaction]
		try:
			rxn = model.reactions.get_by_id(reaction)
			rxn.upper_bound = upperBound
		except KeyError: 
			print reaction

	model.optimize()
	print model.solution.f


	solution = model.solution.x_dict
	ppck = solution['PPCK']
	eno = solution['ENO']
	pepFromOaa = ppck/(ppck+eno)
	lines[0] += str(PEPfromOAA) + "\t"

	tkt2 = solution['TKT2']
	tala = solution['TALA']
	E4PThroughTK = tkt2/(tkt2+tala)
	lines[1] += str(E4PThroughTK) + "\t"

	ppc = solution['PPC']
	mdh = solution['MDH']
	mdh2 = solution['MDH2']
	mdh3 = solution['MDH3']
	OAAfromPEP = ppc/(ppc+mdh+mdh2+mdh3)
	lines[2] += str(OAAfromPEP) + "\t"

	tkt1 = solution['TKT1']
	eda = solution['EDA']
	pfk = solution['PFK']
	PEPthroughTK = tkt1/(eda + pfk + tala + tkt1)
	lines[3] += str(PEPthroughTK) + "\t"

	ghmt2r = solution['GHMT2r']
	serat = solution['SERAT']
	lserdhr = solution['LSERDHr']
	psp_l = solution['PSP_L']
	sert4pp = solution['SERt4pp']
	sert2rpp = solution['SERt2rpp']
	SERfromGLY = ghmt2r/(ghmt2r+serat+lserdhr+psp_l+sert4pp+sert2rpp)
	lines[4] += str(SERfromGLY) + "\t"

	glyat = solution['GLYAT']
	pragsr = solution['PRAGSr']
	amptasecg = solution['AMPTASECG']
	amptasepg = solution['AMPTASEPG']
	sarcox = solution['SARCOX']
	thra2i = solution['THRA2i']
	thrai = solution['THRAi']
	GLYfromSER = ghmt2r/(ghmt2r+glyat+pragsr+amptasecg+amptasepg+sarcox+thra2i+thrai)
	lines[5] += str(GLYfromSER) + "\t"

	try:
		me1 = solution['ME1']
		me2 = solution['ME2']
		pyk = solution['PYK']
		PYRfromMAL = (me1+me2)/(me1+me2+pyk+eda)
		lines[6] += str(PYRfromMAL) + "\t"
	except ZeroDivisionError:
		lines[6] += "0\t"

outFile = open('ConstrainedFluxes_GMPS2_PRFGS_ADSL2r.txt', 'w')
for line in lines:
	outFile.write(line+"\n")
outFile.close()

# wt_solution = model.solution.x_dict
# wt_map = Builder("iJO1366_central_metabolism", reaction_data=wt_solution)
# wt_map.display_in_browser()

	
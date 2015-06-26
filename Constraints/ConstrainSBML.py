"""
	The last script in the sequence. Adds the bounds specified
	in the pickle files in the folder /Vmax. Output two text
	files which have the flux ratios for the specified reactions
	over all the time points, one for the unconstrained fluxes 
	one for the constrained fluxes. 
"""

import pickle
from cobra.io import read_sbml_model
from cobra.core import Model as md
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import os

def changeDefaultBounds (model):

	for rxn in model.reactions:
		if "EX_" in str (rxn):
			if rxn.upper_bound > 1000:
				rxn.upper_bound = 1000
			if rxn.lower_bound < -1000:
				rxn.lower_bound = -1000
	return model

# Equations for flux ratios.
def getFluxRatios (solution):

	solution = solution.x_dict

	try:
		ppck = solution['PPCK']
		eno = solution['ENO']
		PEPfromOAA = ppck/(ppck+eno)
	except ZeroDivisionError:
		PEPfromOAA = 0

	try:
		tkt2 = solution['TKT2']
		tala = solution['TALA']
		E4PThroughTK = tkt2/(tkt2+tala)
	except ZeroDivisionError:
		E4PThroughTK = 0

	try:
		ppc = solution['PPC']
		mdh = solution['MDH']
		mdh2 = solution['MDH2']
		mdh3 = solution['MDH3']
		OAAfromPEP = ppc/(ppc+mdh+mdh2+mdh3)
	except ZeroDivisionError:
		OAAfromPEP = 0 

	try:
		tkt1 = solution['TKT1']
		eda = solution['EDA']
		pfk = solution['PFK']
		PEPthroughTK = tkt1/(eda + pfk + tala + tkt1)
	except ZeroDivisionError:
		PEPthroughTK = 0 

	try:
		ghmt2r = solution['GHMT2r']
		serat = solution['SERAT']
		lserdhr = solution['LSERDHr']
		psp_l = solution['PSP_L']
		sert4pp = solution['SERt4pp']
		sert2rpp = solution['SERt2rpp']
		SERfromGLY = ghmt2r/(ghmt2r+serat+lserdhr+psp_l+sert4pp+sert2rpp)
	except ZeroDivisionError:
		SERfromGLY = 0

	try:
		glyat = solution['GLYAT']
		pragsr = solution['PRAGSr']
		amptasecg = solution['AMPTASECG']
		amptasepg = solution['AMPTASEPG']
		sarcox = solution['SARCOX']
		thra2i = solution['THRA2i']
		thrai = solution['THRAi']
		GLYfromSER = ghmt2r/(ghmt2r+glyat+pragsr+amptasecg+amptasepg+sarcox+thra2i+thrai)
	except ZeroDivisionError:
		GLYfromSER = 0 

	try:
		me1 = solution['ME1']
		me2 = solution['ME2']
		pyk = solution['PYK']
		PYRfromMAL = (me1+me2)/(me1+me2+pyk+eda)
	except ZeroDivisionError:
		PYRfromMAL = 0

	return {'PEPfromOAA': PEPfromOAA, 'E4PThroughTK': E4PThroughTK, 'OAAfromPEP': \
			OAAfromPEP, 'PEPthroughTK': PEPthroughTK, 'SERfromGLY': SERfromGLY,   \
			'GLYfromSER': GLYfromSER, 'PYRfromMAL': PYRfromMAL, 'Growth Rate': model.solution.f}

def constrainModel (model, percent, size, info):

	fluxRatiosList = []

	for time in range (1,10):
		print (percent, size, time)

		fileName = ('Vmax/{}Confidence/{}/Size{}/Time{}.pickle').format(percent, info, size, time)
		inFile = open (fileName, 'rb')
		vMax = pickle.load(inFile)
		inFile.close()

		print ("Changing bounds.")

		count = 0
		for reaction in vMax:
			try:
				rxn = model.reactions.get_by_id(reaction)
				rxn.upper_bound = vMax[reaction]
			except KeyError:
				count += 1

		print ('Could not load {} reactions.').format(count)

		solution = optimize_minimal_flux(model)
		fluxRatios = getFluxRatios (solution)
		fluxRatiosList.append (fluxRatios)
		
	return (fluxRatiosList)

def writeFluxRatios (fluxRatiosList, fileName):

	split = fileName.split('/')
	folderName = ('/').join(split[:-1])
	fileName = ('{}/{}').format(folderName,split[-1])

	# Check if the folder for the file exists, otherwise make it. 
	if not os.path.isdir(folderName):
		os.makedirs (folderName)

	# Make a file in with the file name and add a header. 
	outFile = open (fileName, 'w')
	header = ('\t').join(fluxRatiosList[0]) + '\n'
	outFile.write(header)

	for timepoint in fluxRatiosList:
		line = '\t'.join (str(timepoint[reaction]) for reaction in timepoint) + '\n'
		outFile.write (line)

	outFile.close()

def main ():

	# Read and build SBML REL model. 
	model = read_sbml_model ("REL.xml")
	model = changeDefaultBounds (model)
	solution = optimize_minimal_flux(model)

	# Print growth rate. 
	print (('Growth Rate: {}').format(solution.f))

	# Get the normal flux ratios of the unconstrained model. 
	defaultFluxRatios = getFluxRatios (solution)
	fluxRatiosList = []
	for i in range (9):
		fluxRatiosList.append(defaultFluxRatios)

	writeFluxRatios (fluxRatiosList, 'Results/RegularFluxes.txt')

	# Iterate through each confidence level.
	confidence = [90]
	# Iterate through different list sizes.
	limit = [50]
	# Iterate through each timepoint.
	timepoint = range (9)

	for percent in confidence:
		for size in limit:

			copyModel = md.copy(model)

			fluxRatiosList = constrainModel (copyModel, percent, size, 'RSD')
			fileName = ('Results/{}Confidence/RSD/Size{}.txt').format(percent, size)
			writeFluxRatios (fluxRatiosList, fileName)

			copyModel2 = md.copy(model)

			fluxRatiosList = constrainModel (copyModel2, percent, size, 'Mean')
			fileName = ('Results/{}Confidence/Mean/Size{}.txt').format(percent, size)
			writeFluxRatios (fluxRatiosList, fileName)

main ()


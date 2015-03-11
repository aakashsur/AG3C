"""
	Check growth rate under each condition. 
"""
import cPickle
from cobra.io import *
from cobra.flux_analysis import *
import copy 
# from escher import Builder

def lower_bounds(model):
	print('Lowering bounds.')

	model.optimize()
	newBounds = {999999.000000: 1000.0, -999999.0 : -1000.0}

	for rxn in model.reactions:
		if "EX_" in str(rxn):
			thisReaction = str(rxn)
			try:
				rxn.upper_bound = newBounds[rxn.upper_bound]
				rxn.lower_bound = newBounds[rxn.lower_bound]
			except KeyError:
				pass

	return model

def add_vmax(vMax, model, stop1 = '', stop2 = '', stop3 = ''):

	newModel = copy.deepcopy(model)

	print ("Changing bounds.")

	for reaction in vMax:
		if reaction == stop1 or reaction == stop2 or reaction == stop3:
			continue
		upperBound = vMax[reaction]
		try:
			rxn = newModel.reactions.get_by_id(reaction)
			rxn.upper_bound = upperBound
		except KeyError: 
			pass
			# print (reaction)

	return newModel

def main():

	# Normal Flux
	f1 = open('GrowthRate/RegularFlux.txt', 'w')
	model = read_sbml_model("REL.xml")
	lowerModel = lower_bounds(model)
	lowerModel.optimize()
	solution = lowerModel.solution.f
	print ("Growth Rate is {}").format(solution)
	f1.write(str(solution))
	f1.close()

	vMax = []
	for timepoint in range(0,9):

		# Tells script where to pull pickle.
		fileName = 'Vmax/VMaxTimepoint' + str(timepoint+1) + '.pickle'
		f = open(fileName, "rb")
		vMax.append(cPickle.load(f))
	
	# Calculating fluxes for the constrained models.
	f2 = open('GrowthRate/ConstrainedFlux.txt', 'w')
	f2.write("Timepoint" + '\t' + "Solution" + '\n')

	timepoint = 1

	for bounds in vMax:
		
		model = add_vmax(bounds, lowerModel)
		model.optimize()
		solution = model.solution.f
		print ("Growth Rate at timepoint {} is {}").format(timepoint, solution)
		f2.write(str(timepoint) + '\t' + str(solution) + '\n')
		timepoint += 1

	f2.close()

	f3 = open('GrowthRate/ConstrainedFlux-GMPS2.txt', 'w')
	f3.write("Timepoint" + '\t' + "Solution" + '\n')

	timepoint = 1
	
	for bounds in vMax:

		model = add_vmax(bounds, lowerModel, 'GMPS2')
		model.optimize()
		solution = model.solution.f
		print ("Growth Rate at timepoint {} is {}").format(timepoint, solution)
		f3.write(str(timepoint) + '\t' + str(solution) + '\n')
		timepoint += 1

	f3.close()

	f4 = open('GrowthRate/ConstrainedFlux-GMPS2-PRFGS.txt', 'w')
	f4.write("Timepoint" + '\t' + "Solution" + '\n')

	timepoint = 1
	
	for bounds in vMax:

		model = add_vmax(bounds, lowerModel, 'GMPS2', 'PRFGS')
		model.optimize()
		solution = model.solution.f
		print ("Growth Rate at timepoint {} is {}").format(timepoint, solution)
		f4.write(str(timepoint) + '\t' + str(solution) + '\n')
		timepoint += 1

	f4.close()

	f5 = open('GrowthRate/ConstrainedFlux-GMPS2-PRFGS-ADSL2r.txt', 'w')
	f5.write("Timepoint" + '\t' + "Solution" + '\n')

	timepoint = 1
	
	for bounds in vMax:

		model = add_vmax(bounds, lowerModel, 'GMPS2', 'PRFGS', 'ADSL2r')
		model.optimize()
		solution = model.solution.f
		print ("Growth Rate at timepoint {} is {}").format(timepoint, solution)
		f5.write(str(timepoint) + '\t' + str(solution) + '\n')
		timepoint += 1

	f5.close()

main()

# model = read_sbml_model("REL.xml")
# model = lower_bounds(model)
# model.optimize()
# wt_solution = model.solution.x_dict
# wt_map = Builder("iJO1366_central_metabolism", reaction_data=wt_solution)
# wt_map.display_in_browser()

	
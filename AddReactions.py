"""
This script is designed to modify an existing SBML file to
add, remove, and turn off certain reactions.

Requires a tab delimited file in the format:
"Reaction Abbrevation"  "Reaction"  "Gene"  "Action"

Requires a tab delimited file for details for adding reactions:
"Metabolite Abbrevation"    "Metabolite Name"   "Formula"

Here, I use the iAF1260 SBML model for E. coli MG1655 to develop
a model for E. coli REL.
"""

from cobra.io import *
from cobra.flux_analysis import *
import re

model = read_sbml_model("MG1655.xml")

# Calculate Wild Type Flux Network
model.optimize()
print "Growth Rate: %.15f"%(model.solution.f)

# tab-delimited list of reactions to add, remove, and turn off
f = open("RELReactionOutput.txt", "rU")
read = f.readlines()
f.close()

addReactions = {}
removeReactions = []
turnOffReactions = []

# read the list for the last column which indicates what needs
# to be done with the reaction

for line in read:
    line = line.rstrip()
    split = line.split("\t")

    if split[3] == "deletion":
        removeReactions.append(split[0])

    elif split[3] == "turned off":
        turnOffReactions.append(split[0])

    elif split[3] == "addition": 
        addReactions[split[0]] = split[1]

print '%i reactions in original model' % len(model.reactions)
        
print "Removing reactions."
    
for rxn in removeReactions:
    model.remove_reactions(rxn)

print "Done removing reactions."

print '%i reaction in updated model' % len(model.reactions)

model.optimize()
print "Growth Rate: %.15f"%(model.solution.f)

print 'Turning off reactions.'

for rxn in turnOffReactions:
    try:
        rxn = model.reactions.get_by_id(rxn)
        rxn.upper_bound = rxn.lower_bound = 0
    except:
        rxnCompartment = rxn[-2:-1]
        rxn = rxn[:-3]+"_LPAREN_"+rxnCompartment+"_RPAREN_"
        rxn = model.reactions.get_by_id(rxn)
        rxn.upper_bound = rxn.lower_bound = 0
        
print "Done turning off reactions."

model.optimize()
print "Growth Rate: %.15f"%(model.solution.f)


print "Adding reactions."
print '%i reactions in original model' % len(model.reactions)

# this is the most complex part of parsing

f = open("addReactionsOutput.txt", "rU")
additionalInfo = f.readlines()
f.close()

from cobra import Reaction
from cobra import Metabolite

for rxn in addReactions:

    reaction = Reaction(rxn)
    # Here the reaction name is just the symbolic form of the reaction.
    reaction.name = addReactions[rxn]

    # If the reaction is reversible, set these bounds:
    if "<==>" in addReactions[rxn]:
        lowerBound = -1000
        upperBound = 1000
        
    # If the reaction is directional, set these bounds:
    elif "--" in addReactions[rxn]:
        lowerBound = 0
        upperBound = 1000
    
    reaction.lower_bound = lowerBound
    reaction.upper_bound = upperBound

    # Parse the reactants and products into metabolites.

    rxn = rxn.strip()
    split = re.split(r"(<==>|-->)", addReactions[rxn])
    reactants = split[0].strip()
    products = split[2].strip()

    # If the reaction is an exchange reaction:
    if rxn[:3] == "EX_":
        compartment = "_" + reactants[1:2]
        reactants = reactants[6:].split(" + ")
        products = [reactants[0]]
        reactants[0] += compartment
        products[0] += "_b"

    # If the reaction occurs in one cellular compartment:
    elif reactants[:1] is "[":
        compartment = "_" + reactants[1:2]
        reactants = reactants[6:].split(" + ")
        products = products.split(" + ")
        x = 0
        while x < len(reactants):
            reactants[x] += compartment
            x+=1
        x = 0
        while x < len(products):
            products[x] += compartment
            x+=1

    # If the reaction occurs across multiple compartments:
    else:
        reactants = reactants.split(" + ")
        x = 0
        while x < len(reactants):
            reactants[x] = reactants[x][:-3] + "_" + reactants[x][-2:-1]
            x+=1
        products = products.split(" + ")
        x = 0
        while x < len(products):
            products[x] = products[x][:-3] + "_" + products[x][-2:-1]
            x+=1
        
    """For every reaction object, we must specify the reactant and product
       metabolites. The easist way to do this is if the metabolite is already
       in the model. The script will try to call the metabolite, and will catch
       the error if it does not exist in the model. Next, this program references
       a list of metabolites and their corresponding information. If it still
       cannot find the metabolite, it will throw an error, but not a fatal one.
    """
       
    for molecule in reactants:

        if "-" in molecule:
            molecule = molecule[:-4]+"_DASH_"+molecule[-3:-2]+"_"+molecule[-1:]
        
    
        try:
            reaction.add_metabolites({model.metabolites.get_by_id(molecule): -1})
        except KeyError:
            #print "Error getting reactant %s"%(molecule)
            try:
                mol = None
                for line in additionalInfo:
                    line = line.strip()
                    line = line.split("\t")
                    if molecule[:-2] == line[0].strip():
                        mol = Metabolite(molecule, formula = line[2].strip(), name = line[1].strip(), compartment = molecule[-1:])
                        break
                reaction.add_metabolites({mol:-1})

            except:
                print "Please input data for reactant %s" % (molecule)
            
    for molecule in products:

        if "-" in molecule:
            molecule = molecule[:-4]+"_DASH_"+molecule[-3:-2]+"_"+molecule[-1:]

        try:
            reaction.add_metabolites({model.metabolites.get_by_id(molecule): 1})
        except KeyError:
            #print "Error getting product %s"%(molecule)
            try:
                mol = None
                for line in additionalInfo:
                    line = line.strip()
                    line = line.split("\t")
                    if molecule[:-2] == line[0].strip():
                        mol = Metabolite(molecule, formula = line[2].strip(), name = line[1].strip(), compartment = molecule[-1:])
                        #print mol
                        break
                reaction.add_metabolites({mol:1})

            except:
                print "Please input data for product %s" % (molecule)

    # Final step, add all the reaction object to the SBML model.                
    model.add_reaction(reaction)
    
print '%i reaction in updated model' % len(model.reactions)

model.optimize()
print "Growth Rate: %.15f"%(model.solution.f)

print "writing"
write_sbml_model(model, "REL1.xml")

model = read_sbml_model("REL1.xml")



import cobra.io
import os
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
import numpy as np
# import copy
# from cobra import Metabolite, Reaction
# #import cameo
# import pandas as pd

# functions
def getNamesMediaMetabolites(model):
    listMetabolites = []
    medium = model.medium
    for i in medium:
        mediumReaction = model.reactions.get_by_id(i)
        reactant = mediumReaction.reactants
        nameMetabolite = reactant[0].name
        listMetabolites.append(nameMetabolite)

    return listMetabolites

def exchangeReactions2metaboliteNames(reactionList):
    listMetabolites = []
    reactionDict = {}
    for i in reactionList:
        #mediumReaction = model.reactions.get_by_id(i)
        reactant = i.reactants
        nameMetabolite = reactant[0].name
        listMetabolites.append(nameMetabolite)
        reactionDict.update({nameMetabolite:i})
    return listMetabolites, reactionDict

def productFluxes(model, tol = 1e-6):
    exchangeRnxNames = model.exchanges
    exchangeDict = exchangeReactions2metaboliteNames(exchangeRnxNames)[1]
    solution = model.optimize()
    #tol = 1e-6 # what would be a good threshold? 1e-6 maybe? (tutorial matlab pFBA)
    fluxDict = {}
    for i in exchangeDict: # loop over the name of the
        reactionMetabolite = exchangeDict[i]
        FBA_metabolite_flux = solution.fluxes[reactionMetabolite.id]

        if FBA_metabolite_flux > tol : # threshold for wat is considered as "produced"
            fluxDict.update({i:FBA_metabolite_flux})

    return fluxDict


def plot_flux_solutions(modelLocation, substrate_exchange_rnx, product_exchange_rnx, objectiveReaction = None,
                        conversionFactor = 1.0 , yLabel = 'glucose yield (g/g)' ,FBA = True,  pFBA = True, FVA = True):
    allYields_pFBA =[]
    allYields_FBA =[]
    allYields_FVA_upper = []
    allYields_FVA_lower = []
    modelNames = []
    objectiveBiomass = []
    for i in modelLocation:
        model = cobra.io.read_sbml_model(i)
        # make sure the right objective is set
        if objectiveReaction:
            model.objective = objectiveReaction
        # change the glucose reaction to zero
        glucose_exchange_rnx = 'Ex_S_cpd00027_ext'
        model.reactions.get_by_id(glucose_exchange_rnx).bounds = 0,0
        # change bound of new substrate to -10 mol/h/gDW
        model.reactions.get_by_id(substrate_exchange_rnx).bounds = -10, 0
        # get names of the models
        modelName = i.split("\\")[-1]
        modelName = modelName.replace(".xml","")
        modelNames.append(modelName)
        # run pFBA
        if pFBA:
            pfba_solution = cobra.flux_analysis.pfba(model)
            substrate_flux = pfba_solution.fluxes[substrate_exchange_rnx]
            product_flux = pfba_solution.fluxes[product_exchange_rnx]
            yield_pFBA = -product_flux/substrate_flux # fixed
            allYields_pFBA.append(yield_pFBA)
            # get objective flux
            #objectiveBiomass.append(pfba_solution.objective_value) #this gives a wierd result
            objectiveBiomass.append(pfba_solution.fluxes["Ex_S_biomass_ext"])

        if FVA:
            bounds_fva = flux_variability_analysis(model,reaction_list= [substrate_exchange_rnx,product_exchange_rnx])
            bounds_yield_lower = -bounds_fva['minimum'][product_exchange_rnx]/bounds_fva['minimum'][substrate_exchange_rnx]
            bounds_yield_upper = -bounds_fva['maximum'][product_exchange_rnx]/bounds_fva['maximum'][substrate_exchange_rnx]
            allYields_FVA_lower.append(bounds_yield_lower)
            allYields_FVA_upper.append(bounds_yield_upper)

        if FBA:
            solution = model.optimize()
            FBA_substrate_flux = solution.fluxes[substrate_exchange_rnx]
            FBA_product_flux = solution.fluxes[product_exchange_rnx]
            FBA_yield = -FBA_product_flux/FBA_substrate_flux
            allYields_FBA.append(FBA_yield)
            objectiveBiomass.append(solution.objective_value)

    # conversionFactor =  92.1/180 propionate/glycerol
    # conversionFactor = 74/180 # propionate/glucose to go from mol to grams change to one for moles
    bottom_ = np.multiply(allYields_FVA_lower, conversionFactor)
    height_ = np.multiply(allYields_FVA_upper, conversionFactor)
    x_cordinates = range(len(modelLocation))
    plt.bar(x = x_cordinates, height=height_, width=0.8, bottom=bottom_,tick_label = modelNames)
    if pFBA:
        plt.plot(x_cordinates,np.multiply(allYields_pFBA,conversionFactor),'y*', markersize = 18)
    if FBA:
        fbas = np.multiply(allYields_FBA, conversionFactor)
        print(fbas)
        plt.plot(x_cordinates, fbas,'r*', markersize = 10)
    plt.grid()
    plt.ylabel(yLabel)
    plt.legend(['pFBA','FBA','FVA'])
    plt.ylim(bottom= 0)

    return plt, objectiveBiomass


if __name__ == '__main__':

    loc = os.getcwd()
    loc_acidi = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML models\PAC_4875_model.xml' #loc + r'\SBML models\PAC_4875_model.xml'
    loc_acnes = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_acnes_model.xml'
    loc_prop = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_propionicum_model.xml'
    loc_avidum = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_avidum_model.xml'
    loc_sher = r'C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_sherm_model.xml'

    #microorganisms = [loc_acidi, loc_acnes, loc_prop]
    microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher]
    #microorganisms = [loc_acnes, loc_sher]

    glucose_exchange_rnx = 'Ex_S_cpd00027_ext'
    propionate_exchange_rnx = 'Ex_S_cpd00141_ext'
    acetate_exchange_rnx = 'Ex_S_cpd00029_ext'
    fructose_exchange_rnx = 'Ex_S_cpd00082_ext'
    glycerol_exchange_rnx = 'Ex_S_cpd00100_ext'
    # model.exchanges[0].id # is the objective reaction for propioni models


    plotFig = True
    if plotFig:
        pltGlu, objetiveGlu = plot_flux_solutions(microorganisms, substrate_exchange_rnx=glucose_exchange_rnx,
                                                  product_exchange_rnx=propionate_exchange_rnx, conversionFactor=3/6,
                                                  pFBA=False, yLabel= 'glucose-prop yields ')
        print('glucose-prop yields ')
        #print(objetiveGlu)
        pltGlu.show()

        pltGlu, objetiveGlu = plot_flux_solutions(microorganisms, substrate_exchange_rnx=glucose_exchange_rnx,
                                                  product_exchange_rnx=acetate_exchange_rnx, conversionFactor=2/6,
                                                  pFBA=False, yLabel= 'glucose-ace yields ')
        print('glucose-ace yields ')
        # print(objetiveGlu)
        pltGlu.show()

        pltGlu, objetiveGlu = plot_flux_solutions(microorganisms, substrate_exchange_rnx=fructose_exchange_rnx,
                            product_exchange_rnx=propionate_exchange_rnx,conversionFactor= 3/6, pFBA= False,
                                                  yLabel= 'frutose-prop yields ')
        print('fructose-prop yields ')
        # print(objetiveGlu)
        pltGlu.show()

        pltGlu, objetiveGlu = plot_flux_solutions(microorganisms, substrate_exchange_rnx=fructose_exchange_rnx,
                            product_exchange_rnx=acetate_exchange_rnx,conversionFactor= 2/6, pFBA= False,
                                                  yLabel= 'frutose-acetate yields ')

        print('fructose-acetate yields ')
        # print(objetiveGlu)
        pltGlu.show()

        # pltGlu, objetiveGlu = plot_flux_solutions(microorganisms, substrate_exchange_rnx=glucose_exchange_rnx,
        #                     product_exchange_rnx=propionate_exchange_rnx, objectiveReaction= "Ex_S_biomass_ext"
        #                     ,conversionFactor= 74/180)
        #print(objetiveGlu)
        #pltGlu.show()

        # # # # # # # # #  glycerol plot# # # # # # # # # # # #
        # pltGly, objetiveGly = plot_flux_solutions(microorganisms,substrate_exchange_rnx=glycerol_exchange_rnx,
        #                                 product_exchange_rnx=propionate_exchange_rnx, conversionFactor = 92.1/180,
        #                                 objectiveReaction= "Ex_S_biomass_ext", yLabel= 'yield glycerol (g/g)')
        # print(objetiveGly)
        # pltGly.show()
         # # # # # # # # #  glycerol plot# # # # # # # # # # # #

    testArea = False
    if testArea:  #handy functions and ways to call metabolites or reactions
        # use acicdi species as an example
        model = cobra.io.read_sbml_model(loc_acidi)
        test = model.exchanges
         # look for metabolites in media
        mediumNames = getNamesMediaMetabolites(model=model) # names in medium
        print(mediumNames)

        # look for the sinks, demand or exchange/boundry reactions:
        # more info: https: // cobrapy.readthedocs.io / en / latest / media.html
        exchangeRnxNames = model.exchanges # An exchange reaction is a reversible reaction that adds to or removes an # extracellular metabolite from the extracellular compartment.
        boundryRnxNames = model.boundary
        sinkRnxNames = model.sinks # A sink is similar to an exchange but specifically for intracellular metabolites, i.e., a reversible reaction that adds or removes an intracellular metabolite.
        demandRnxNames = model.demands # A demand reaction is an irreversible reaction that consumes an intracellular metabolite.

        # Find name that are interpretable
        exchangeMet = exchangeReactions2metaboliteNames(exchangeRnxNames)[0]
        boundryMet = exchangeReactions2metaboliteNames(boundryRnxNames)[0]
        sinkMet = exchangeReactions2metaboliteNames(sinkRnxNames)[0]
        demandMet = exchangeReactions2metaboliteNames(demandRnxNames)[0]

        # print results
        print(len(exchangeMet))
        print(exchangeMet)
        print("")
        print(len(boundryMet))
        print(boundryMet)
        print("")
        print(sinkMet)
        print(demandMet)

        # want to know if an exchange/boundry reaction is consumed or produced? -> check bounds and run a FBA
        #for example, wat happens to Methanol?
        exchangeDict = exchangeReactions2metaboliteNames(exchangeRnxNames)[1]
        reactionMethanol = exchangeDict['Methanol']
        print(reactionMethanol.bounds)
        #bound : [0, 1000] so does not get consumed?
        # so does it get produced? run FBA
        solution = model.optimize()
        FBA_methanol_flux = solution.fluxes[reactionMethanol.id]
        print(FBA_methanol_flux) # no flux, so we can treat this metabolite as a substrate by changeing the bounds!!
        # check out the functions to see how to get names with  model.reaction.get_from_id()

        # or use slef made function
        # finds produced products of the exchange reactions
        model = cobra.io.read_sbml_model(loc_acidi)
        fluxes = productFluxes(model)
        print(fluxes)

        ## look for specifice reaction/ metabolite
        # for example you want to know wat Ex_S_cpd00001_ext is (water by the way):
        model.reactions.get_by_id('Ex_S_cpd00001_ext') #and look into object while in debug mode



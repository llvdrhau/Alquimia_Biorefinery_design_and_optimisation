
import cobra.io
import os
from cobra.flux_analysis import flux_variability_analysis
import matplotlib.pyplot as plt
# import copy
# from cobra import Metabolite, Reaction
# #import cameo
# import pandas as pd

# functions

def plot_flux_solutions(modelLocation, substrate_exchange_rnx, product_exchange_rnx,FBA = True, pFBA = True, FVA = True):
    allYields_pFBA =[]
    allYields_FBA =[]
    allYields_FVA_upper = []
    allYields_FVA_lower = []
    modelNames = []
    for i in modelLocation:
        model = cobra.io.read_sbml_model(i)
        model.reactions.get_by_id(substrate_exchange_rnx).bounds = -10,0
        modelName = i.split("\\")[-1]
        modelName = modelName.replace(".xml","")
        modelNames.append(modelName)
        # run pFBA
        if pFBA:
            pfba_solution = cobra.flux_analysis.pfba(model)
            substrate_flux = pfba_solution.fluxes[substrate_exchange_rnx]
            product_flux = pfba_solution.fluxes[product_exchange_rnx]
            yield_pFBA =  -substrate_flux/product_flux
            allYields_pFBA.append(yield_pFBA)

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

    bottom_ = allYields_FVA_lower
    height_ = allYields_FVA_upper
    x_cordinates = range(len(modelLocation))
    plt.bar(x = x_cordinates, height=height_, width=0.8, bottom=bottom_,tick_label = modelNames)
    plt.plot(1,2,'b*')
    plt.show()

if __name__ == '__main__':

    loc = os.getcwd()
    loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
    loc_acnes = loc + r'\SBML models\P_acnes_model.xml'
    loc_prop = loc + r'\SBML models\P_propionicum_model.xml'

    #microorganisms = [loc_acidi, loc_acnes, loc_prop]
    microorganisms = [loc_acnes, loc_prop]

    glucose_exchange_rnx = 'Ex_S_cpd00027_ext'
    propionate_exchange_rnx = 'Ex_S_cpd00141_ext'
    plot_flux_solutions(microorganisms,substrate_exchange_rnx=glucose_exchange_rnx,
                       product_exchange_rnx=propionate_exchange_rnx)

    testArea = False
    if testArea:
        #model = cobra.io.read_sbml_model('../SBML models/PAC_4875_model.xml')
        model = cobra.io.read_sbml_model(loc_acidi)
        #
        # # reactions of interest
        # biomass_rnx = 'biomass_c0'
        # glucose_exchange_rnx = 'Ex_S_cpd00027_ext'
        # propionate_exchange_rnx = 'Ex_S_cpd00141_ext'
        # # uptake rate of glucose bound to 10 mol/h/g_dw
        # model.reactions.Ex_S_cpd00027_ext.boounds = (-10,0)
        # # run pFBA
        # pfba_solution = cobra.flux_analysis.pfba(model)
        # glucose_flux = pfba_solution.fluxes[glucose_exchange_rnx]
        # propionate_flux = pfba_solution.fluxes[propionate_exchange_rnx]
        #
        # bounds = flux_variability_analysis(model, propionate_exchange_rnx)
        #
        # cv = propionate_flux/glucose_flux
        #
        # print(model.exchanges)
        # print(cv)
        # print(bounds)
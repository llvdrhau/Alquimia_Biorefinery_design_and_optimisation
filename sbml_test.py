
import cobra.io
import os
# import copy
from cobra.flux_analysis import flux_variability_analysis
# from cobra import Metabolite, Reaction
# #import cameo
# import pandas as pd

loc = os.getcwd()
loc = loc + r'\SBML models\PAC_4875_model.xml'

# model = cobra.io.read_sbml_model('../SBML models/PAC_4875_model.xml')
model = cobra.io.read_sbml_model(loc)

# reactions of interest
biomass_rnx = 'biomass_c0'
glucose_exchange_rnx = 'Ex_S_cpd00027_ext'
propionate_exchange_rnx = 'Ex_S_cpd00141_ext'
# uptake rate of glucose bound to 10 mol/h/g_dw
model.reactions.Ex_S_cpd00027_ext.boounds = (-10,0)
# run pFBA
pfba_solution = cobra.flux_analysis.pfba(model)
glucose_flux = pfba_solution.fluxes[glucose_exchange_rnx]
propionate_flux = pfba_solution.fluxes[propionate_exchange_rnx]

bounds = flux_variability_analysis(model, propionate_exchange_rnx)

cv = propionate_flux/glucose_flux

print(model.exchanges)
print(cv)
print(bounds)
"""
fixing the all propionibacteria models in the same manner

lucas.vanderhauwaert@usc.es
18/07/2023
"""

from function_fix_P_models import fix_P_models_sbml

# from the script C:\Users\Lucas\PycharmProjects\Alquimia\SBML screening\Fixing_models\atp_maintanance_Bacillus.py
# we can see that the maintenance ATP is 9 mmol/gDW/h for B. subtilis which is a gram positive bacteria like propionibacteria
# and a facultative anaerobe like propionibacteria (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC106906/)
# so we will use the same value for the maintenance ATP for all propionibacteria models

fix_P_models_sbml(modelName='P_avidum_model.xml',      boundATPm= 9, saveName='P_avidum_V2.xml')
fix_P_models_sbml(modelName='P_propionicum_model.xml', boundATPm= 9, saveName='P_propionicum_V2.xml')
fix_P_models_sbml(modelName='P_acnes_model.xml',       boundATPm= 9, saveName='P_acnes_V2.xml')
fix_P_models_sbml(modelName='P_sherm_model.xml',       boundATPm= 9, saveName='P_sherm_V2.xml')
fix_P_models_sbml(modelName = 'PAC_4875_model.xml',    boundATPm= 9, saveName='PAC_4875_V2.xml')
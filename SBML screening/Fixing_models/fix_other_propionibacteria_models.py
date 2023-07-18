"""
fixing the all propionibacteria models in the same manner

lucas.vanderhauwaert@usc.es
18/07/2023
"""

from function_fix_P_models import fix_P_models_sbml

fix_P_models_sbml(modelName='P_avidum_model.xml', boundATPm= 3.15, saveName='P_avidum_V2.xml')
fix_P_models_sbml(modelName='P_propionicum_model.xml', boundATPm= 3.15, saveName='P_propionicum_V2.xml')
fix_P_models_sbml(modelName='P_acnes_model.xml', boundATPm= 3.15, saveName='P_acnes_V2.xml')
fix_P_models_sbml(modelName='P_sherm_model.xml', boundATPm= 3.15, saveName='P_sherm_V2.xml')
fix_P_models_sbml(modelName = 'PAC_4875_model.xml', boundATPm= 3.15, saveName='PAC_4875_V2.xml')

from f_usefull_functions import get_location
from f_screen_SBML import find_yield, string_reactions
import cobra

import cobra
from mewpy.simulation import get_simulator

if __name__ == '__main__':
    # list of models to be analyzed
    ListModels = [ 'P_acnes_V2.xml', 'PAC_4875_V2.xml', 'P_propionicum_V2.xml','P_avidum_V2.xml',  'P_sherm_V2.xml']
    modelName = ListModels[-1]

    locationModel = get_location(modelName)
    model = cobra.io.read_sbml_model(locationModel)

    # reaction IDs for to examin the FVB results
    reactions = [
        #propionate
        'Ex_S_cpd00141_ext',
        #acetate
        'Ex_S_cpd00029_ext', ]

    simul = get_simulator(model)  # replace with your model's file path
    fva_results = simul.FVA(reactions=reactions, obj_frac=0.95)
    print(fva_results)
    print('')

    propionate_yield_min = fva_results['Ex_S_cpd00141_ext'][0] * 73 / (10* 180.16)    # 10 mmol/gDW/h of substrate
    propionate_yield_max = fva_results['Ex_S_cpd00141_ext'][1] * 73 / (10* 180.16)      # 10 mmol/gDW/h of substrate
    acetate_yield_max = fva_results['Ex_S_cpd00029_ext'][1]* 60.05     / (10* 180.16)           # 10 mmol/gDW/h of substrate
    acetate_yield_min = fva_results['Ex_S_cpd00029_ext'][0]* 60.05     / (10* 180.16)                 # 10 mmol/gDW/h of substrate

    # get the yield
    glucoseRnxId = 'Ex_S_cpd00027_ext'
    propionate_yield_fba = find_yield(model=model, substrateExchangeRxnID=glucoseRnxId,
                                      productExchangeRxnID='Ex_S_cpd00141_ext', optimisationMode='FBA')
    propionate_yield_pfba = find_yield(model=model, substrateExchangeRxnID=glucoseRnxId,
                                      productExchangeRxnID='Ex_S_cpd00141_ext', optimisationMode='pFBA')
    # print the results
    print('The model is: ', modelName)
    print('')
    print('The minimum propionate yield is: ', propionate_yield_min)
    print('The maximum propionate yield is: ', propionate_yield_max)
    print('The FBA propionate yield is: {}'.format(propionate_yield_fba))
    print('The pFBA propionate yield is: {} \n'.format(propionate_yield_pfba))

    print('The minimum acetate yield is: ', acetate_yield_min)
    print('The maximum acetate yield is: ', acetate_yield_max)




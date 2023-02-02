import cobra.io
import os
import re


def find_carbons_in_formula(formula):
    metFormula = formula
    splitFormula = re.split('(\d+)', metFormula)
    nrOfCarbons = 0  # just in case something wierd happens
    if 'C' not in metFormula:  # if there is no carbon in the formula
        nrOfCarbons = 0
    else:
        for j, element in enumerate(splitFormula):
            if 'C' in element and len(element) == 1:
                nrOfCarbons = int(splitFormula[j + 1])
            elif 'C' in element and len(element) > 1:
                posCarbon = element.index('C')
                if element[posCarbon + 1].isupper():  # for case like CH4 there is no 1 next to the C
                    nrOfCarbons = 1  # only one carbon
                else:
                    continue  # for cases like Co (cobalt) just skip
            else:
                continue
    return nrOfCarbons

def get_conversion_sbml(modelLocations, substrate_exchange_rnx, product_exchange_rnx, substrate2zero= 'Ex_S_cpd00027_ext',
                        newObjectiveReaction = None, pFBA = None, printEq = False):
    allYields_pFBA =[]
    allYields_FBA =[]
    objectiveBiomass = []
    allEquations = []
    modelNames = []
    for i in modelLocations:
        model = cobra.io.read_sbml_model(i)
        # make sure the right objective is set
        if newObjectiveReaction:
            model.objective = newObjectiveReaction
        # change the glucose reaction to zero
        exchange_rnx_2_zero = substrate2zero
        model.reactions.get_by_id(exchange_rnx_2_zero).bounds = 0,0

        # change bound of new substrate to -10 mol/h/gDW
        model.reactions.get_by_id(substrate_exchange_rnx).bounds = -10, 0
        # get names of the models
        modelName = i.split("\\")[-1]
        modelName = modelName.replace(".xml","")
        modelNames.append(modelName)
        # run pFBA
        if pFBA: # Todo fix flux to grams c for pFBA
            pfba_solution = cobra.flux_analysis.pfba(model)
            substrate_flux = pfba_solution.fluxes[substrate_exchange_rnx]
            product_flux = pfba_solution.fluxes[product_exchange_rnx]
            yield_pFBA = -product_flux/substrate_flux # fixed
            allYields_pFBA.append(yield_pFBA)

        else: #else do regular FBA
            solution = model.optimize()
            FBA_substrate_flux = solution.fluxes[substrate_exchange_rnx]
            substrateMet = model.reactions.get_by_id(substrate_exchange_rnx).reactants[0]
            substrateName = substrateMet.name
            substrateFormula = substrateMet.formula
            Csub = find_carbons_in_formula(substrateFormula)
            strEqlist = []

            for j in product_exchange_rnx:
                productMet = model.reactions.get_by_id(j).reactants[0]
                productName = productMet.name
                productFormula = productMet.formula
                Cprod = find_carbons_in_formula(productFormula)
                FBA_product_flux = solution.fluxes[j]
                FBA_yield = abs((FBA_product_flux/FBA_substrate_flux) * (Cprod *12) /(Csub* 12)) # in gramsC / grams C: 12 gCarbon/mol

                allYields_FBA.append(FBA_yield)
                strEq = '{} == {} * {}'.format(productName,FBA_yield,substrateName)
                allEquations.append(strEq)
                if printEq:
                    print(modelName)
                    print(strEq)

    return allEquations, allYields_FBA

if __name__ == '__main__':
    loc = os.getcwd()
    loc_acidi = loc + r'\SBML models\PAC_4875_model.xml'
    loc_acnes = loc + r'\SBML models\P_acnes_model.xml'
    loc_prop = loc + r'\SBML models\P_propionicum_model.xml'
    loc_avidum = loc + r'\SBML models\P_avidum_model.xml'
    loc_sher = loc + r'\SBML models\P_sherm_model.xml'
    microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum, loc_sher] # all microorganisms

    products = ['Ex_S_cpd00029_ext', 'Ex_S_cpd00141_ext'] # acetate and propionate
    substrates = ['Ex_S_cpd00027_ext', 'Ex_S_cpd00082_ext']
    #aa = get_conversion_sbml(modelLocations= microorganisms,substrate_exchange_rnx= 'Ex_S_cpd00027_ext', product_exchange_rnx= products, printEq= True)
    #a = get_conversion_sbml2(modelLocation= loc_acidi , substrate_exchange_rnx= substrates, product_exchange_rnx= products, printEq=True, checkCarbon=True)
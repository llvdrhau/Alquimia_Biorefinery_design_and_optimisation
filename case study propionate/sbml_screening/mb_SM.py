
import cobra
import numpy as np
from f_find_carbons import countCarbonInFormula

# todo try and make sens of the stoichiometrix matrix got dmn it

def find_metabolite_index(model, metID):
    mets = model.metabolites
    for i, met in enumerate(mets):
        if met.id == metID:
            break
    return i

loc_sher = r"C:\Users\lucas\PycharmProjects\Alquimia\SBML models\P_sherm_model.xml"
model = cobra.io.read_sbml_model(loc_sher)
stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(model)
FBA = model.optimize()
fluxArray = FBA.fluxes.to_numpy()

metabolicFlux = np.dot(stoichiometric_matrix, fluxArray)
metabolite_carbon_flux_array = []
positive = []
negative = []
for i,met in enumerate(model.metabolites):
    nCarbon = countCarbonInFormula(met.formula)
    carbonFlux = nCarbon * metabolicFlux[i]* 12
    metabolite_carbon_flux_array.append(carbonFlux)
    if carbonFlux > 0:
        positive.append(carbonFlux)
    elif carbonFlux < 0:
        negative.append(carbonFlux)

gluMetIndex = find_metabolite_index(model,'S_cpd00027_ext')
glu_row = stoichiometric_matrix[gluMetIndex,:]

print(sum(metabolite_carbon_flux_array))
print('the ratio in out is: {}'.format(abs(sum(positive)/sum(negative))))

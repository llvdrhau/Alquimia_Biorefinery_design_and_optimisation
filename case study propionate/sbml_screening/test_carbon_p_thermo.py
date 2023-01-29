from print_model_2_excel import print_SBML_info_2_excel
from f_usefull_functions import get_location
import cobra


modelName = "p-thermo.xml"
loc_sher = get_location(modelName)
model = cobra.io.read_sbml_model(loc_sher)
model.name = modelName

print_SBML_info_2_excel(modelName=model)
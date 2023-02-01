
from old_scripts.old2.f_find_carbons import carbonBalanceInOut
import cobra
import os
cdw = os.getcwd()
print(cdw)


loc = os.getcwd()
microViviane = loc + r'\SBML models\p-thermo.xml'  #iNF517.xml' #
validat = cobra.io.sbml.validate_sbml_model(microViviane)
print(validat[0])
print('')
print(validat[1]['SBML_FATAL']) #
#model = cobra.io.read_sbml_model(microViviane)
carbonBalanceInOut(modelLocation=microViviane, metIDsMissingCarbon=[], tol=1e-4)


from f_screen_SBML import find_yield, string_reactions
from f_usefull_functions import get_location
import cobra

modelName = "iYO844.xml"
modelLocation = get_location(modelName)
model = cobra.io.read_sbml_model(modelLocation)


ATPMid="ATPM"
atpmRxn = model.reactions.get_by_id(ATPMid)
strATPm = string_reactions(atpmRxn, printFlux=True)
print('')
print("The maintenance reaction is:")
print(strATPm[0])
print('the bounds of the reaction are: {}'.format(atpmRxn.bounds))
print('The flux of the maintanance reaction is: {}'.format(strATPm[1]))



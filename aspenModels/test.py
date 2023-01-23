

import os
import win32com.client as win32
# aspen = win32.Dispatch('Apwn.Document')
aspen = win32.gencache.EnsureDispatch("Apwn.Document")
modelPath = os.path.abspath(r'aspenModels\test_model.bkp')
aspen.InitFromArchive2(modelPath)

tree_path = r"\Data\Blocks\B2\Input\NSTAGE"
trays = aspen.Tree.FindNode(tree_path).Value
print('the amount of trays are:')
print(trays)

aspen.Close()

# #access to the specific node
# node = aspen.Tree.FindNode(tree_path)
# #get the measure unit of the node (codified)
# node_unit = node.AttributeValue(ATTRNAME_MAP['Unit'])
# #get the node basis (Mole,Mass,…)
# node_basis = node.AttributeValue(ATTRNAME_MAP['Basis'])
# node.SetValueUnitAndBasis(Value=value, unitcol=node_unit, basis=node_basis)
#



# from Code_Library import Simulation
#
# sim = Simulation(AspenFileName= "test_model.bkp",
#                  WorkingDirectoryPath= r"C:\Users\lucas\PycharmProjects\Alquimia\aspenModels"
#                  ,VISIBILITY=False)
# sim.CloseAspen()
# print('model opended and closed')
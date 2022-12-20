import cobra
import numpy as np
import pandas as pd

test_model = cobra.Model("test_model")
v1 = cobra.Reaction("v1")
v2 = cobra.Reaction("v2")
v3 = cobra.Reaction("v3")
v4 = cobra.Reaction("v4")
v5 = cobra.Reaction("v5")
v6 = cobra.Reaction("v6")
vextra = cobra.Reaction("vextra")

test_model.add_reactions([v1, v2, v3, v4, v5, v6, vextra])

v1.reaction = "-> 2 A"
v2.reaction = "A <-> B"
v3.reaction = "A -> D"
v4.reaction = "A -> C"
v5.reaction = "C -> D"
v6.reaction = "D ->"
vextra.reaction = "B -> x + 11 C"

v1.bounds = (-10, 10)
v2.bounds = (-10.0, 10.0)
v3.bounds = (0.0, 3.0)
v4.bounds = (0.0, 3.0)
v5.bounds = (0.0, 3.0)
v6.bounds = (0.0, 3.0)
vextra.bounds = (0.0, 10.0)

test_model.objective = v6

FBA = test_model.optimize()
fluxArray = FBA.fluxes #[0:posExchangeRxn] #.to_numpy()
print(fluxArray)

stoichiometric_matrix = cobra.util.create_stoichiometric_matrix(test_model, array_type='DataFrame')  # [:,0:posExchangeRxn]
print(stoichiometric_matrix)
# stoichiometric_matrix.drop(ExRxn, axis=1, inplace=True)  # get rid of the exchange reactions
# stoichiometric_matrix.drop(ExMet, axis=0, inplace=True)  # get rid of the exchanged metabolites

metabolicFlux = np.dot(stoichiometric_matrix, fluxArray)
metabolicFluxDF = pd.DataFrame(metabolicFlux, index=stoichiometric_matrix.index)
print(metabolicFluxDF)
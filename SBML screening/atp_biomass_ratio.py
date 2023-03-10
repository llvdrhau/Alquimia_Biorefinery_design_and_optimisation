from f_screen_SBML import *

modelList = ["P_sherm_model.xml", 'P_propionicum_model.xml', 'PAC_4875_model.xml', 'P_acnes_model.xml', 'P_avidum_model.xml']

for modelName in modelList:
    # read the model
    loc_sher = get_location(modelName)
    model = cobra.io.read_sbml_model(loc_sher)

    # get the flux solution
    solutions = model.optimize()
    fluxArray = solutions.fluxes

    # get the stoichiometric matrix
    StoiMatrixDF = cobra.util.create_stoichiometric_matrix(model, array_type='DataFrame')  # [:,0:posExchangeRxn]

    # metabolite and reaction id's
    biomassExReactionID = 'Ex_S_biomass_ext'
    ATPMetaboliteID = 'S_cpd00002_c0'

    # find the ATP producing flux
    RowATP = StoiMatrixDF.loc[ATPMetaboliteID]
    fluxATP = RowATP * fluxArray

    fluxList = []
    for flux in fluxATP:
        if flux > 0:
            fluxList.append(flux)
    totalATPFlux = sum(fluxList)

    # find the biomass flux
    biomassExReaction = model.reactions.get_by_id(biomassExReactionID)
    biomassFlux = biomassExReaction.flux

    # the ratio
    ATP_BM_ratio = totalATPFlux / biomassFlux

    print('model:', modelName)
    print('the total flux of ATP is: ', totalATPFlux)
    print('the total flux of BM is: ', biomassFlux)
    print('the total ratio ATP/BM is: ', ATP_BM_ratio)
    print('')





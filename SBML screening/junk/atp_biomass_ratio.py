from f_screen_SBML import *

modelList = ["P_sherm_lvdh.xml", 'P_propionicum_lvdh.xml', 'PAC_4875_lvdh.xml', 'P_acnes_lvdh.xml', 'P_avidum_lvdh.xml']

for modelName in modelList:
    # read the model
    loc_sher = get_location(modelName)
    model = cobra.io.read_sbml_model(loc_sher)

    # make aerobic
    # O2RxnID = 'Ex_S_cpd00007_ext'
    # O2Rnx = model.reactions.get_by_id(O2RxnID)
    # O2Rnx.bounds = -10, 0

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
    biomassFlux = biomassExReaction.flux        # in mmol/gDW/h
    metBM = biomassExReaction.reactants[0]
    formula_BM = metBM.formula
    MW_BM =  metBM.formula_weight               # in mg/mmol (equivalent to g/mol)
    massFluxBiomass = biomassFlux * MW_BM # in mg/gDW/h

    # the ratio
    BM_ATP_ratio = massFluxBiomass / totalATPFlux  # in mg/mmol or g/mol

    print('model:', modelName)
    print('the MW of biomass is: {} \n the formula of biomass is: {} '.format(MW_BM,formula_BM ))
    print('the total flux (mmol/gDW/h) of ATP is: ', totalATPFlux)
    print('the total flux (mg/gDW/h) of BM is: ', biomassFlux)
    print('the total ratio BM/ATP (g/mol) is: ', BM_ATP_ratio)
    print('')





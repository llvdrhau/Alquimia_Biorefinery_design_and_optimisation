import pandas as pd
from f_make_surrogate_model import *

regresionSurrogate = True
SBMLsurrogate = False

# -------------------------- open fermentation
if regresionSurrogate:
    # Glucose_PH_Data.xlsx # without biomass as variable
    # excelFile = 'Glucose_PH_Data.xlsx'
    excelFile = 'Glucose_PH_60_Data_Points.xlsx'  # with a 5th polynomial it is acctually quite nice....

    # ---- plot the data
    dataLocation = get_location(file=excelFile, case='ML')
    x_pH = pd.read_excel(dataLocation, sheet_name= 'inputs')
    y_outputs = pd.read_excel(dataLocation, sheet_name= 'outputs')
    #plot_subplots(y_data=y_outputs, x_data=x_pH)

    # deleet the row where the pH is larger than 8.49
    # (otherwise the fit is not going to be great, see the propionate plot)

    # get the indexes of rows with values greater than 8.48
    deleted_indexes = x_pH[x_pH < 8.5].index
    x_pH = x_pH[x_pH < 8.5]
    y_outputs.drop(deleted_indexes)
    plot_subplots(y_data=y_outputs, x_data=x_pH)


    # ---- fit poly data
    polynomial = 4
    out = regression_2_json(excelFile, normalise= False ,save=False, saveName='open_fermentation_polynomial_case_study.json',
                         showPLot= True, polynomial= {'pH':polynomial}, case= 'Ridge')

# -------------------------- SBML models
if SBMLsurrogate:

    loc_acidi = 'PAC_4875_lvdh.xml'
    loc_acnes = 'P_acnes_lvdh.xml'
    loc_prop = 'P_propionicum_lvdh.xml'
    loc_avidum = 'P_avidum_lvdh.xml'
    loc_sher = 'P_sherm_lvdh.xml'

    # SBML surrogates
    substrates = ['Ex_S_cpd00027_ext', # glucose
                  'Ex_S_cpd00082_ext', # fructose
                   'Ex_S_cpd00100_ext' # glycerol
                  ]

    products = ['Ex_S_cpd00141_ext',   # acetate
                'Ex_S_cpd00029_ext',    # propionate
                'Ex_S_biomass_ext']     # biomass

    microorganisms = [loc_acidi, loc_acnes, loc_prop, loc_avidum]  # all microorganisms
    saveNames = ['PAC.json', 'acnes.json', 'propionicum.json' ,'avidum.json']
    for i, organism in enumerate(microorganisms):
        sbml = SBML_2_json(modelName= organism, substrate_exchange_rnx= substrates, product_exchange_rnx= products,
                           saveName= saveNames[i], save = False, checkCarbon= False, printEq= True)

    # the model P_sherm_model can not consume glycerol so only screen it for glucose and fructose 'sherm'

    organism = loc_sher
    saveName = 'sherm.json'

    substrates = ['Ex_S_cpd00027_ext', # glucose
                  'Ex_S_cpd00082_ext', # fructose
                  ]

    products = ['Ex_S_cpd00141_ext',   # acetate
                'Ex_S_cpd00029_ext',    # propionate
                'Ex_S_biomass_ext']     # biomass

    sbmlSherm = SBML_2_json(modelName=organism, substrate_exchange_rnx=substrates, product_exchange_rnx=products,
                           saveName=saveName, save=False, checkCarbon=False, printEq=True)
"""
this scripts contains the functions to make the reaction equations as strings
for the surrogate models
"""

from f_get_conversion_from_SBML import get_conversion_sbml
import os
# ============================================================================================================
# Functions to make the interval reaction equations from SBML models
# ============================================================================================================
def split_remove_spaces(expr2split,splitCharacter):
    expresions = expr2split.split(splitCharacter)
    exprList = []
    for i in expresions:
        i = i.replace(' ', '')
        exprList.append(i)
    return exprList

def make_str_eq_smbl(modelName, substrate_exchange_rnx, product_exchange_rnx, equationInfo):
    #modelLocations = modelName
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    #modelLocations = loc + r'\excel files\' + modelName
    modelLocations = [loc + r"\SBML models\{}".format(modelName)]

    equations, yields = get_conversion_sbml(modelLocations, substrate_exchange_rnx, product_exchange_rnx)

    # make abbreviations dictionary
    abbrDict = {}
    inputSbmlName = equationInfo.input_name
    inputAbrr = equationInfo.input_abrr
    abbrDict.update({inputSbmlName:inputAbrr})

    outputSbmlName = split_remove_spaces(equationInfo.output_name,',')
    outputAbrr = split_remove_spaces(equationInfo.output_abrr, ',')
    for i, out in enumerate(outputSbmlName):
        abbrDict.update({out:outputAbrr[i]})

    allEquations = []
    for eq in equations:
        for name in abbrDict:
            if name in eq:
                eq = eq.replace(name, abbrDict[name])
        allEquations.append(eq)

    return allEquations
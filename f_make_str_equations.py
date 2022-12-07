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

    # TODO make extended reaction equtions for sbml models or maybe even make them into JSON files? idk
    # would run quicker in stead of having to read the xlm files each time
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

def make_str_eq_json(modelObject, equationInfo):
    outputs = modelObject['outputs']
    coef = modelObject['coef']
    intercept = modelObject['intercept']
    name = modelObject['name']
    lable = modelObject['lable']
    yieldOf = equationInfo['yield_of']
    # make abbreviations dictionary
    abbrDict = {}
    varName = split_remove_spaces(equationInfo.input_name, ',') + split_remove_spaces(equationInfo.output_name, ',')
    Abrr = split_remove_spaces(equationInfo.input_abrr, ',') + split_remove_spaces(equationInfo.output_abrr, ',')
    for i, varN in enumerate(varName):
        abbrDict.update({varN: Abrr[i]})

    equationList = []
    print(name)
    for out in outputs:
        outAbrr = abbrDict[out]
        #eqY = '{} == '.format(outAbrr)
        coefOfOutputs = coef[out]
        yieldEq = ''
        for feature in coefOfOutputs:
            ### this for loop to replace all the full variable names with the abbreviuations
            featureAbbr = ''
            for v in varName:
                if v in feature:
                    featureAbbr = feature.replace(v,abbrDict[v])
            if not featureAbbr:
                raise Exception('the feature {} has no abbreviation check the JSON file and the sheet models, '
                                'are names and abrr correct?'.format(feature))
            ###
            featureCoef = coefOfOutputs[feature]
            yieldEq += ' + {} * {} '.format(featureAbbr,featureCoef)
        if lable == 'SBML':
            equation = '{} == {}'.format(outAbrr, yieldEq)
        else:
            yieldEq += ' + {}'.format(intercept[out])
            equation = '{} == ({}) * {}'.format(outAbrr,yieldEq,yieldOf)
        print(equation)
        print('')
        equationList.append(equation)

    return equationList



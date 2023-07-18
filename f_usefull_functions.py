import os
import json

########################################################################################################################
# ============================================================================================================
# Usefull functions
# ============================================================================================================
########################################################################################################################

def split_remove_spaces(expr2split, splitCharacter):
    """ Splits a given string according to a specified character e.g., ','
    also removes the spaces ' '
    """
    exprList = []
    if not isinstance(expr2split, str) and not isinstance(expr2split, list):  # so a float or int just one number
        return [float(expr2split)]  # if mulptiple prices for utilty you're gona have to change str to floats

    elif isinstance(expr2split, list):
        for exp in expr2split:
            expresions = exp.split(splitCharacter)
            for i in expresions:
                i = i.replace(' ', '')  # remove annoying spaces
                exprList.append(i)
        return exprList

    elif isinstance(expr2split, str):
        expresions = expr2split.split(splitCharacter)
        for i in expresions:
            i = i.replace(' ', '')  # remove annoying spaces
            exprList.append(i)
        return exprList


def stringbounds_2_tuplebounds(stringBound):
    """ transforms a bound that is writen as a string (e.g., '[20, 50]') to a tuple
    """
    stringBound = stringBound.replace('[', '')
    stringBound = stringBound.replace(']', '')
    bounds = stringBound.split(',')
    boundsList = []
    for i in bounds:
        boundsList.append(float(i))
    return boundsList


def load_objectes_from_dictionary(dict):
    """ deleet I think, not used"""
    for i in dict:
        locals()[i] = dict[i]


def remove_spaces(listOfInterest):
    """ removes spaces """
    exprList = []
    for i in listOfInterest:
        i = i.replace(' ', '')
        exprList.append(i)
    return exprList


def get_connected_intervals(intervalName, conectionMatrix):
    """From the conection matrix (see Excel file) the connected intervals are found which go into the specified
    interval (intervalName) is found
    """
    conectionCol = conectionMatrix[intervalName]
    connectionInfo = conectionCol.where(conectionCol != 0).dropna()

    # drop the bool variable (otherwise mixing gets confused)
    # if there is a boolean varibale at least
    try:
        connectionInfo = connectionInfo.drop(intervalName)
    except:
        pass

    nameConnectedIntervals = connectionInfo.index
    connectionDict = {nameConnectedIntervals[i]: connectionInfo[i] for i in range(len(connectionInfo))}

    # posConnect = conectionCol != 0
    # nameConnectedIntervals = list(conectionMatrix['process_intervals'][posConnect])

    return connectionDict


def str_2_dict(string, intervalname):
    D = eval(string)
    inputBoundsDict = {}
    for i in D:
        inputBoundsDict.update({i + '_' + intervalname: D[i]})
    return inputBoundsDict


def get_location(file, case = ''):
    """ gets the file location from the Directory 'Excel files'
    """
    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]


    if '/' in loc:  # in the case of macOS
        file = r"/{}".format(file)
        if '.xlsx' in file:
            if case == 'ML': # if you want to get Exel files for machinelearing go to other location
                loc = loc + r'/machine learning models/excel_data' + file
            else:
                loc = loc + r'/excel files' + file

        elif '.xml' in file:
            loc = loc + r'/SBML models' + file
        elif '.json' in file:
            loc = loc + r'/json models' + file


    elif "\\" in loc:  # in the case of Windows OS
        file = r"\{}".format(file)
        if '.xlsx' in file:
            if case == 'ML': # if you want to get Exel files for machinelearing go to other location
                loc = loc + r'\machine learning models\excel_data' + file
            else:
                loc = loc + r'\excel files' + file

        elif '.xml' in file:
            loc = loc + r'\SBML models' + file
        elif '.json' in file:
            loc = loc + r'\json models' + file

    return loc

def save_2_json(saveName, saveObject):
    """ saves an object or dictionary to a json file """

    if not isinstance(saveObject, dict):
        saveObject = saveObject.__dict__

    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia + 8]
    loc = loc + r'\json models' + r'\{}'.format(saveName)
    with open(loc, 'w+', encoding='utf-8') as f:
        json.dump(saveObject, f, ensure_ascii=False, indent=4)


def linear_aprox(y1, y2, x1, x2, z1):
    v = x1 - (x1-x2)*(y1-z1)/(y1-y2)
    return v


def transform_dictionary(input_dict):
    output_dict = {}

    for key1, value1 in input_dict.items():
        for key2, value2 in value1.items():
            if key2 not in output_dict:
                output_dict[key2] = {}
            output_dict[key2][key1] = value2

    return output_dict

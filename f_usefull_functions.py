# ============================================================================================================
# Usefull functions
# ============================================================================================================
def split_remove_spaces(expr2split,splitCharacter):
    exprList = []
    if not isinstance(expr2split, str) and not isinstance(expr2split,list): # so a float or int just one number
        return [float(expr2split)]   # if mulptiple prices for utilty you're gona have to change str to floats

    elif isinstance(expr2split,list):
        for exp in expr2split:
            expresions = exp.split(splitCharacter)
            for i in expresions:
                i = i.replace(' ', '')  # remove annoying spaces
                exprList.append(i)
        return exprList

    elif isinstance(expr2split, str):
        expresions = expr2split.split(splitCharacter)
        for i in expresions:
            i = i.replace(' ', '') # remove annoying spaces
            exprList.append(i)
        return exprList

def stringbounds_2_tuplebounds(stringBound):
    stringBound = stringBound.replace('[', '')
    stringBound = stringBound.replace(']', '')
    bounds = stringBound.split(',')
    boundsList = []
    for i in bounds:
        boundsList.append(float(i))
    return boundsList

def load_objectes_from_dictionary(dict):
    for i in dict:
        locals()[i] = dict[i]

def remove_spaces(listOfInterest):
    exprList = []
    for i in listOfInterest:
        i = i.replace(' ','')
        exprList.append(i)
    return exprList

def get_connected_intervals(intervalName,conectionMatrix):
    conectionCol = conectionMatrix[intervalName]
    posConnect = conectionCol != 0
    nameConnectedIntervals = list(conectionMatrix['process_intervals'][posConnect])
    connectionInfo = list(conectionMatrix[intervalName][posConnect])
    connectionDict = {nameConnectedIntervals[i]:connectionInfo[i] for i in range(len(connectionInfo))}
    return  connectionDict

def str_2_dict(string,intervalname):
    D = eval(string)
    inputBoundsDict = {}
    for i in D:
        inputBoundsDict.update({i+'_'+intervalname : D[i]})
    return inputBoundsDict
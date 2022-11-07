'''
Created on Tue Oct 04 2022
Contains the classes to make the process interval objects

@author: Lucas Van der Hauwaert
email: lucas.vanderhauwaert@usc.es
'''

# ============================================================================================================
# Input, reactor and output Classes
# ============================================================================================================

class InputIntervalClass:
    def __init__(self, inputName, compositionDict, inputPrice, boundry, boolDict = None, split=None, separationDict=None):
        if separationDict is None:
            separationDict = {}
        if split is None:
            split = []
        if boolDict is None:
            boolDict = {}

        # declare input interval name
        self.label = 'input'
        self.inputName = inputName #.upper() # nah don't put in capitals
        self.inputPrice = inputPrice
        addOn4Variables = '_' + inputName.lower()

        # change the composition names in the dictionary
        compositionDictNew = {}
        initialCompositionNames = []
        for i in compositionDict:
            compositionDictNew.update({i + addOn4Variables: compositionDict[i]})
            initialCompositionNames.append(i)
        self.initialCompositionNames = initialCompositionNames
        self.compositionDict = compositionDictNew

        # error if you choose a wrong name
        for i in initialCompositionNames:  # maybe don't place the warning here
            if i in inputName:
                raise Exception("the component {} is in interval name {}, change the component name of the reactor to "
                                "avoid conflict with the equations")


        # make the component equations as string equations
        eqList = []
        for component in compositionDictNew:
            eq = "{} == {} * {}".format(component, self.compositionDict[component], self.inputName)
            eqList.append(eq)
        self.componentEquations = eqList

        componentVariables =  list(compositionDictNew.keys()) # [self.inputName] not gona add this atm
        # make boundary for all variables
        boundaryDict = {self.inputName:boundry} # intervals have bounds
        for i in componentVariables:
            boundaryDict.update({i:(0, None)})
        self.boundaries = boundaryDict
        self.boolDict = boolDict
        self.split = split
        self.separationDict = separationDict
        #self.compositionNames = list(compositionDict.keys())
        self.leavingInterval = componentVariables

        eqSumOfBoolsHelp = '1 == '
        eqSumOfBools = []
        if boolDict:
            for interval in boolDict:
                eqSumOfBoolsHelp += ' + ' + boolDict[interval]
            eqSumOfBools = [eqSumOfBoolsHelp]
        self.eqSumOfBools = eqSumOfBools

    #  put all VARIABLES that pyomo needs to declare here
        self.componentVariables = componentVariables
        self.boolVariables = list(boolDict.values())

        self.allVariables = componentVariables + [self.inputName] + self.boolVariables

    #  put all EQUATIONS that pyomo needs to declare here
    # self.allEquations =

    def makeOldNewDict(self, oldNames, newNames):
        oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
        self.oldNewDictOutputs = oldNewDict

    def UpdateDict(self, newKeys):
        oldDict = self.compositionDict
        oldKeys = list(oldDict.keys())
        rangeKeys = range(len(oldKeys))
        for i in rangeKeys:
            new_key = newKeys[i]
            old_key = oldKeys[i]
            self.compositionDict[new_key] = self.compositionDict.pop(old_key)

class ReactorIntervalClass:
    def __init__(self, inputs, outputs, reactionEquations, nameDict, mix=None, utilities=None, boolActivation = None,
                 boolDict=None, split=None, separationDict=None):
        if utilities is None:
            utilities = {}
        if separationDict is None:
            separationDict = {}
        if split is None:
            split = []
        if boolDict is None:
            boolDict = []
        if mix is None:
            mix = []
        if boolActivation is None:
            boolActivation = []

        self.label = 'reactor'
        self.nameDict = nameDict
        self.inputs = inputs
        self.outputs = outputs
        self.initialCompositionNames = inputs + outputs
        # error if you choose a wrong name
        for i in self.initialCompositionNames: # maybe don't place the warning here
            if i in list(nameDict.keys())[0]:
                raise Exception("the component {} is in interval name {}, change the component name of the reactor to "
                                "avoid conflict with the equations")

        self.mix = mix  # found by the Excel file (don't need to specify in this script)
        self.utilities = utilities # consists of a dictionary {nameUtilty: [bounds]}
        self.boolDict = boolDict  # is a tuple, 1) where the bool comes from and 2) the name of the bool affecting the outputs
        self.split = split  # TODO will probs be a dict
        self.separationDict = separationDict  # dictionary defining separation fractions of each component

        # interval Variable name
        intervalVariable = list(self.nameDict.keys())[0]

        # reactor equations
        if isinstance(reactionEquations, str):  # failsafe if you forget to code the string reactor expression as a list
            reactionEquations = [reactionEquations]  # make it a list
        else:
            pass
        ouputs2change  = self.outputs
        allEquations = []
        reactionVariablesOutput = []
        rctVarOutD = {}
        for eq in reactionEquations:
            for out in ouputs2change:
                newOutputName = out + '_{}'.format(intervalVariable)
                eq = eq.replace(out,newOutputName)
                if newOutputName not in reactionVariablesOutput:
                    reactionVariablesOutput.append(newOutputName)
                    rctVarOutD.update({out:newOutputName}) # helpìng dictionary for other
            allEquations.append(eq)
        self.reactionEquations = allEquations

        # mass equations (of the outputs from the reaction equations )
        eqMassInterval  = intervalVariable + " == "
        for out in reactionVariablesOutput:
            eqMassInterval += " + " + out
        self.totalMassEquation =  [eqMassInterval]

        # bool activation constraints
        boolActivationEquations = []
        if boolActivation: # if there is an activation constraint
            bounds = nameDict[intervalVariable]
            lowerActivationEq = "{} * {} <= {}".format(boolActivation[0],bounds[0],intervalVariable)
            upperActivationEq = "{} <= {} * {}".format(intervalVariable,boolActivation[0], bounds[1])
            boolActivationEquations.append(lowerActivationEq)
            boolActivationEquations.append(upperActivationEq)
        self.boolActivationEquations =  boolActivationEquations

        # bounds of the component and interval variable
        boundaryDict = {}  # intervals have bounds
        if not boolActivationEquations: # if the reactor is not affected by the boolean constraint you can add it to the boundry list
            boundInterval = nameDict[intervalVariable]
            boundaryDict.update({intervalVariable:boundInterval})
        for i in reactionVariablesOutput:
            boundaryDict.update({i: (0, None)})
        self.boundaries = boundaryDict

        # sum of bool equations
        eqSumOfBoolsHelp = '1 == '
        eqSumOfBools = []
        if self.boolDict:
            for interval in boolDict:
                eqSumOfBoolsHelp += ' + ' + boolDict[interval]
            eqSumOfBools = [eqSumOfBoolsHelp]
        self.eqSumOfBools = eqSumOfBools

        # separation equations
        separationEquations = []
        separationVaribles = []
        for sep in separationDict: # if it is empty it should not loop nmrly
            for componentSep in separationDict[sep]:
                var = rctVarOutD[componentSep]  # helpìng dictionary to get the right variable
                sepVar  = componentSep + '_' + sep
                eqSep = sepVar + " == {} * {} ".format(separationDict[sep][componentSep], var)
                separationEquations.append(eqSep)
                separationVaribles.append(sepVar)

        self.seprationEquations = separationEquations
        self.separationVaribles = separationVaribles

        # spliting equations
        # TODO

        # mixing equations
        # see def make_mixing_equations()

        # define wat is leaving the reactor
        # can either be from the reactor, the separation process or the spliting
        if not split and not separationDict:
            self.leavingInterval = reactionVariablesOutput
        elif separationDict:
            self.leavingInterval = separationVaribles
        elif split:  # todo next!!
            pass
        # todo combination of split and serparation after a reaction, how to do that?

        # put all self.VARIABLES that pyomo needs to declare here
        self.reactionVariablesOutput = reactionVariablesOutput
        # reactionVariablesInputs can be found in class function: update_reactor_equations
        self.intervalVariable = intervalVariable
        if boolActivation: # if it is present
            self.activationVariable = boolActivation[0]
        #make a list whith all of the variables

    def make_mix_equations(self, objects2mix):
        mixEquations = []
        initialInputNames = self.inputs
        intervalNames2Mix = objects2mix.keys()

        # find the leaving variables
        leavingVars = []
        for objName in objects2mix:
            obj = objects2mix[objName][0] # first element in tuple is the object
            connectInfo = objects2mix[objName][1] # first element in tuple is connection info

            # if there is a separation process
            if isinstance(connectInfo, str) and 'sep' in connectInfo:
                allSepVars = obj.leavingInterval
                for var in allSepVars:
                    if connectInfo in var:
                        leavingVars.append(var)
            # otherwise just add the leaving variables
            else:
                leavingVars += obj.leavingInterval

        # make the equations
        mixingVariables = []
        for i, ins in enumerate(initialInputNames):
            mixVar = ins + "_mix"
            eqMix = mixVar + " == "
            startMixEq = eqMix
            for lvar in leavingVars:
                if ins in lvar:
                    eqMix += " + " + lvar
            """
            For example in the case of pH this does not come from the previous reactor!! 
            so the variable can stay as it is and no extra equations needs to be added, hence if eqMix != startMixEq:  
            """
            if eqMix != startMixEq:
                mixEquations.append(eqMix)
                mixingVariables.append(mixVar)
        self.mixEquations = mixEquations
        self.mixingVariables = mixingVariables

        # now the reaction equations need to be updated because
        # they're the mix variables are now the inputs of the reaction equations!
        # see function update_reactor_intervals & class function update_reactor_equations

    def update_reactor_equations(self, newInputs4Interval):
        initialInputs4Interval = self.inputs
        replacementDict = self.get_replacement_dict(initialInputs4Interval, newInputs4Interval)
        equationsInterval = self.reactionEquations  # the reactor equations
        allEquations = []
        for eq in equationsInterval:
            for var in replacementDict:
                newVarName = replacementDict[var]
                eq = eq.replace(var, newVarName)
            allEquations.append(eq)
        self.reactionEquations = allEquations
        self.reactionVariablesInputs = list(replacementDict.values())

    def get_replacement_dict(self,initialVars, newVars):
        replacementDict = {}
        for i in initialVars:
            if i == 'pH':  # TODO find a better way to do this: pH always stays pH_intervalName
                replacementDict.update({i: 'pH_{}'.format(self.intervalVariable)})
            for j in newVars:
                if i in j:  # the initial variable (of excel) is always in the new name, that's how you can find wat belongs to where
                    replacementDict.update({i: j})

        return replacementDict

    # def makeDict(self):
    #     in_outDict = {
    #         'inputs': self.inputs,
    #         'outputs': self.outputs
    #     }
    #     return in_outDict
    #
    # def make_replacement_dict_output(self, oldNames, newNames):
    #     oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
    #     self.replacementDictOutput= oldNewDict
    #     return oldNewDict
    #
    # def make_replacement_dict_input(self, oldNames, newNames):
    #     oldNewDict = {oldNames[i]: newNames[i] for i in range(len(oldNames))}  # make dictionary
    #     self.replacementDictInput = oldNewDict
    #     return oldNewDict
    #
    # def makeSeperationDict(self, permeate,
    #                        reject):  # define dictionary to see wat percentage of each stream component goes to the permeate and reject streams
    #     percentExtractionVector = list(self.separation.values())
    #     seperationDict = {self.outputs[i]: [permeate[i], reject[i], percentExtractionVector[i]] for i in
    #                       range(len(self.outputs))}
    #     return seperationDict
    #
    #
    #     #  could you define a stream dependent on a boolean var before in enters a unit reactor

class OutputIntervalClass:
    def __init__(self, outputName, outputBound ,outputPrice, outputVariable, mixDict = None):
        if mixDict is None:
            mixDict = {}
        self.label = 'output'
        self.outputName = outputName
        self.outputPrice = outputPrice
        self.outputBound = outputBound
        self.outputVariable = outputVariable  # eg: acetate is output name but is refered to as ace in all the reactions
    def make_end_equations(self, objects2mix):

        # find the leaving variables of the object(s) that go into the output
        leavingVars = []
        for objName in objects2mix:
            obj = objects2mix[objName][0] # first element in tuple is the object
            connectInfo = objects2mix[objName][1] # first element in tuple is connection info

            # if there is a separation process
            if isinstance(connectInfo, str) and 'sep' in connectInfo:
                allSepVars = obj.leavingInterval
                for var in allSepVars:
                    if connectInfo in var and self.outputVariable in var:
                        leavingVars.append(var)
            # otherwise just add the leaving variables
            else:
                leavingVars += obj.leavingInterval

        # make the equations
        endVar = self.outputName
        eqEnd = endVar + " == "
        #startMixEq = eqEnd
        for i, lvar in enumerate(leavingVars):
            eqEnd += " + " + lvar

        self.endEquations = eqEnd
        self.endVariables = endVar

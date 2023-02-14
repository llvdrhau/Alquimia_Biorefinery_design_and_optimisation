from functions import *

test = False
if test:
    excelFile = r'\propionate_case_study.xlsx'

    loc = os.getcwd()
    posAlquimia = loc.find('Alquimia')
    loc = loc[0:posAlquimia+8]
    loc = loc + r'\excel files' + excelFile
    DFConnect = pd.read_excel(loc, sheet_name='connection_matrix_test', index_col='process_intervals')

    DFConnect = DFConnect.drop(labels= [ 'carbon_input1', 'carbon_input2',
                                        'carbon_input3', 'acetate','propionate','waste'], axis=1)
    interValNames = [
    'P_acidi',
    'P_freu',
    'P_avi',
    'P_acn',
    'P_pro',
    'open_fermentation',
    'liq_liq_ext',
    'NF',
    'Distilation_1',
    'Distilation_2',
    'Distilation_3']


    booleanVariables, equationsSumOfBools = make_boolean_equations(DFConnect,interValNames)
    print(booleanVariables)
    print(equationsSumOfBools)

# def get_helping_dict_4_seperation(self, connectInfo):
    #     """ makes the helping dictionary to make the seperation equations.
    #     IMPORTANT: the variables are not declared here!! this happens in the update functions. This is a hack to get the
    #     variables you need based on how the interval is connected to the previous intervals.
    #     SEE t_oDO in the function update_intervals
    #     e.g., {'ace': 'ace_P_freu_batch_sep1', 'prop': 'prop_P_freu_batch'_sep1, 'water': 'water_P_freu_batch_sep1'}
    #
    #     Parameters:
    #         connectInfo (dict)
    #
    #     Returns:
    #         helpDict (dict)
    #         """
    #     #  fix this conundrum:
    #     #  the variables for mixing are declared in the update functions (called after this function!!) but here we are
    #     #  replicating those varibles. In other words if you change the mixing, seperation or split variable names this
    #     #  part of the code will certainly give an error... NEED TO FIND A CONSISTENT WAY TO CALL INCOMING STREAMS
    #     #   Idealy get rid of the update functions and make all equations immidiatly
    #
    #     originalOutputNames = self.outputs
    #     intervalName = list(self.nameDict.keys())[0]
    #     helpDict = {}
    #     if len(connectInfo) >1:  # then the components are mixed
    #         # entering vars are then mixed
    #         for outName in originalOutputNames:
    #             enteringVar = "{}_{}_mix".format(outName,intervalName)
    #             helpDict.update({outName:enteringVar})
    #
    #     elif len(connectInfo) == 1: #components are separated
    #         intervalKey = list(connectInfo.keys())[0]
    #         connectSpecification = connectInfo[intervalKey] #
    #         reactorKey ,sepKey, splitKey , boolKey = define_connect_info(connectSpecification)
    #         if splitKey:
    #             for outName in originalOutputNames:
    #                 enteringVar = "{}_{}_{}".format(outName, intervalName,splitKey)
    #                 helpDict.update({outName: enteringVar})
    #
    #         elif sepKey:
    #             for outName in originalOutputNames:
    #                 enteringVar = "{}_{}_{}".format(outName, intervalName,sepKey)
    #                 helpDict.update({outName: enteringVar})
    #
    #         else:
    #             pass # then it is just from the raeactor which the function make_reaction_equations already generates
    #
    #     return helpDict



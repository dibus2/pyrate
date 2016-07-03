#!/usr/bin/env python 
try:
    import pudb
    import yaml
    import sys
    import argparse
    import os
    import time
    wd = os.getcwd()
    sys.path.append(wd + '/src/Core')
except:
    exit("Error while importing one of the modules `sys,os,argparse,yaml`")

welcomemessage = """\n\n\n\n\t\t==================================================================================\n
\t\t\t\tPyR@TE version 2.0.0  released  ???\n
\t\t\tF. Lyonnet, I. Schienbein,\n
\t\t\tand F.Staub, A.Wingerter (version 1)
\t\t==================================================================================\n
"""
print(welcomemessage)
# Read the Model file and create the class
# Create the parser
parser = argparse.ArgumentParser(description='')
# Add an argument model i.e. the name of the YAML file with the settings
parser.add_argument('--Settings', '-f', dest='Settings', action='store', default='', help='Specify a settings file')
parser.add_argument('--Model', '-m', dest='Model', action='store', default='',
                    help="This option is used to specify the YAML file containing the model settings")
parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', default=False,
                    help='Prints details on the calculation')
parser.add_argument('--VerboseLevel', '-vL', dest='VerboseLevel', action='store', default='',
                    help='Set up the level of printing')
parser.add_argument('--Gauge-Couplings', '-gc', dest='Gauge-Couplings', action='store_true', default=False,
                    help='Set the calculation to the gauge couplings beta functions')
parser.add_argument('--Quartic-Couplings', '-qc', dest='Quartic-Couplings', action='store_true', default=False,
                    help='Set the calculation to the quartic couplings beta functions')
parser.add_argument('--Yukawas', '-yu', dest='Yukawas', action='store_true', default=False,
                    help='Set the calculation to the yukawas beta functions')
parser.add_argument('--ScalarMass', '-sm', dest='ScalarMass', action='store_true', default=False,
                    help='Set the calculation to the Scalar Mass beta functions')
parser.add_argument('--FermionMass', '-fm', dest='FermionMass', action='store_true', default=False,
                    help='Set the calculation to the Fermion Mass beta functions')
parser.add_argument('--Trilinear', '-tr', dest='Trilinear', action='store_true', default=False,
                    help='Set the calculation to the Trilinear term beta functions')
parser.add_argument('--All-Contributions', '-a', dest='All-Contributions', action='store_true', default=False,
                    help='Set the calculation to all the contributions')
parser.add_argument('--Two-Loop', '-tl', dest='Two-Loop', action='store_true', default=False,
                    help='Set the calculation to two loop order')
parser.add_argument('--LogFile', '-lg', dest='LogFile', action='store', default='logFile.log',
                    help='Set the logFile name')
parser.add_argument('--LogLevel', '-lv', dest='LogLevel', action='store', default='Info', help='Set the LogLevel ')
parser.add_argument('--LatexFile', '-tex', dest='LatexFile', action='store', default='RGEsOutput.tex',
                    help='Set the name of the Latex output file')
parser.add_argument('--LatexOutput', '-texOut', dest='LatexOutput', action='store_true', default=True,
                    help='Write the results in a Latex file')
parser.add_argument('--Results', '-res', dest='Results', action='store', default='./results/',
                    help='store all the results in the path')
parser.add_argument('--Pickle', '-pkl', dest='Pickle', action='store_true', default=False,
                    help='Set if a pickle output should be produced')
parser.add_argument('--PickleFile', '-pf', dest='PickleFile', action='store', default='RGEsOutput.pickle',
                    help='Set the name of the pickle output file')
parser.add_argument('--TotxtMathematica', '-tm', dest='ToM', action='store_true', default=False,
                    help='Set if a pickle output should be produced')
parser.add_argument('--TotxtMathFile', '-tmf', dest='ToMF', action='store', default='RGEsOutput.txt',
                    help='Set the name of the txt to mathematica file')
parser.add_argument('--Weyl', '-w', dest='Weyl', action='store', default=True,
                    help='Set if the fields are weyl or dirac, i.e. if kappa should be taken equal to 1 or 1/2')
parser.add_argument('--Export', '-e', dest='Export', action='store_true', default=False,
                    help='Set if a python function is created to evaluate the beta function')
parser.add_argument('--Export-file', '-ef', dest='ExportFile', action='store', default='BetaFunction.py',
                    help='Set the name of the output python file')
parser.add_argument('--Only', '-onl', dest='Only', action='store', default=[],
                    help='Set a dictionary of terms you want to calculate: "QuarticTerms,Yukawas,TrilearTerms,FermionMasses,ScalarMasses". E.g. "{\'QuarticTerms\': [\'\lambda_1\',\'\lambda_2\']}". Note that if passed in the command line the whole argument must be a string')
parser.add_argument('--Skip', '-sk', dest='Skip', action='store', default='',
                    help='Set the different terms to neglect in the calculation. E.g. ["CAabcd","CL2abcd"]. The list of terms that can be neglected are listed in src/Core/RGEsDefinition.py')
parser.add_argument('--Pylie', '-ipl', dest='interactivedb', action='store_true', default=False,
                    help='Starts the interactive PyLie mode. Allows one to check what are the CGCs implemented for a given contraction and more (Casimir, Dynkin,...)')
parser.add_argument('--ScalarAnomalous', '-sa', dest='ScalarAnomalous', action='store_true', default=False,
                    help='Set the calculation of scalar anomalous dimensions to True')
parser.add_argument('--FermionAnomalous', '-fa', dest='FermionAnomalous', action='store_true', default=False,
                    help='Set the calculation of fermion anomalous dimensions to True')
parser.add_argument('--SetGutNorm', '-gutn', dest='SetGutNorm', action='store_true', default=False,
                    help='Set the normalization to gut normalization in case there is a U(1) gauge group, it normalizes g1 -> sqrt(3/5)*g\'')
parser.add_argument('--KinMix', '-kin', dest='KinMix', action='store_false', default=True,
                    help='Switch off the kinetic mixing terms if multiple U(1) gauge groups are present.')

# Collect the arguments
args = parser.parse_args()
FileSettings = args.__getattribute__('Settings')
# Default Run Settings on top of the default in the command line
# Copy the argument either from the settings File or from the command line
if FileSettings == "":
    # Transform args into a dictionnary and copy it
    RunSettings = args.__dict__
else:
    # Copy the default values
    RunSettings = {}
    try:
        f = open('{}'.format(FileSettings), 'r')
        yamlRunSettings = yaml.load(f.read())
        f.close()
    except yaml.scanner.ScannerError as err:
        exit("Check the YAML file {}, failed to load the settings: {}".format(FileSettings, err))
    except yaml.parser.ParserError as err:
        exit("Check the YAML file {}, failed to parse the settings: {}".format(FileSettings, err))
    except IOError as (errno, errstr):
        exit(
            "Did not find the YAML file {}, specify the path if not in the current directory.\n {}".format(FileSettings,
                                                                                                           errstr))
    try:
        for key, value in yamlRunSettings.items():
            RunSettings[key] = value
    except KeyError as err:
        exit("EROR, key {} does not exist, refer to the list of defined keys:\n\n {}".format(key, err))
    # overwrite the settings by the commmand line
    defaults = dict([(el, parser.get_default(el)) for el in args.__dict__.keys()])
    for key, val in args.__dict__.items():
        if key in RunSettings and val != defaults[key]:
            RunSettings[key] = val
        elif not (key in RunSettings):
            RunSettings[key] = val
        else:
            pass
try:
    import copy
    sys.path.append(wd + '/src/Output')
except ImportError:
    exit("Error while importing the modules")
from Logging import *
# get current working directory
cwd = os.getcwd()
RunSettings['Results'] = os.path.abspath(RunSettings['Results'])
if not (os.path.exists(RunSettings['Results'])):
    os.makedirs(RunSettings['Results'])
# Get the logLevel from RunSettings
LogLevel = {'Info': logging.INFO, 'Debug': logging.DEBUG, 'Critical': logging.CRITICAL}
# create the config of the logging system
if RunSettings['LogFile'] == 'logFile.log':
    # add the time stamp
    RunSettings['LogFile'] = RunSettings['LogFile'].split('.')[0] + '_{}.log'.format(time.time())

logging.basicConfig(filename=RunSettings['LogFile'], level=LogLevel['{}'.format(RunSettings['LogLevel'])],
                    format="%(levelname)s [%(asctime)s] [%(funcName)s] %(message)s")
# Setting up the verbose level
if RunSettings['verbose']:
    if RunSettings['VerboseLevel'] == 'Info':
        RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = True, False, True
    elif RunSettings['VerboseLevel'] == 'Debug':
        RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = False, True, True
    elif RunSettings['VerboseLevel'] == 'Critical':
        RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = False, False, True
    elif RunSettings['VerboseLevel'] == '':
        RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = True, False, True
    else:
        print("unknown VerboseLevel {}...setting it to Info".format(RunSettings['VerboseLevel']))
        RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = True, False, False
else:
    RunSettings['vInfo'], RunSettings['vDebug'], RunSettings['vCritical'] = False, False, True

# Interrogate the database
if 'interactivedb' in RunSettings and RunSettings['interactivedb']:
    from IPyrate import *
    Idbquerry().cmdloop()
    loggingInfo("\nExiting the interactive mode.", verbose=True)
    exit()

# Set the skipping terms
if 'Skip' in RunSettings:
    if type(RunSettings['Skip']) == str:
        RunSettings['Skip'] = RunSettings['Skip'][1:-1].split(',')
# For the loggingInfo function the RunSettings dict can be passed to print on the screen if verbose is true by default it is true in loggingInfo
loggingInfo("Starting a new run ...", verbose=RunSettings['vInfo'])
strRunSettings = '\n'.join(["\t\t{}: {}".format(key, val) for key, val in RunSettings.items()])
loggingInfo("Run Settings :\n{}".format(strRunSettings), verbose=RunSettings['vInfo'])
# Processing of the required argument
if RunSettings['Model'] == "":
    loggingCritical("Please, specify a .model file.", verbose=RunSettings['vCritical'])
    exit()
else:
    try:
        # open the YAML file
        f = open('{}'.format(RunSettings['Model']), 'r')
        # Open the Yaml file and load the settings
        yamlSettings = yaml.load(f.read())
        f.close()
    except yaml.scanner.ScannerError as err:
        loggingCritical(
            "Check the YAML file {}, impossible to load the settings: {}.".format(RunSettings['Model'], err),
            verbose=RunSettings['vCritical'])
        exit()
    except yaml.parser.ParserError as err:
        loggingCritical(
            "Check the YAML file {}, impossible to parse the settings: {}.".format(RunSettings['Model'], err),
            verbose=RunSettings['vCritical'])
        exit()

    except IOError as (errno, errstr):
        loggingCritical("Did not find the YAML file {}, specify the path if not in the current directory.\n {}.".format(
            RunSettings['Model'], errstr), verbose=RunSettings['vCritical'])
        exit()
    loggingInfo("Loading the YAML file: {} ... done".format(RunSettings['Model']), verbose=RunSettings['vInfo'])

    # Now we want to process the settings before creating the model class
    # Let's first construct the dictionaries if the input is given as a list
    if 'Groups' in yamlSettings and yamlSettings['Groups'] != {} and type(yamlSettings['Groups']) == list:
        tp = {}
        cg = []
        ListG = []
        for x in yamlSettings['Groups']:
            if yamlSettings['Groups'].count(x) == 1:
                tp.update({'{}'.format(x): x})
                ListG.append(x)
            else:
                cg.append(x)
                tp.update({'{}_{}'.format(x, cg.count(x)): x})
                ListG.append('{}_{}'.format(x, cg.count(x)))
        yamlSettings['Groups'] = copy.deepcopy(tp)
    if 'Fermions' in yamlSettings and yamlSettings['Fermions'] != {} and type(
            yamlSettings['Fermions'].values()[0]) == list:
        if not (all([len(pp) == len(yamlSettings['Groups']) + 1 for pp in yamlSettings['Fermions'].values()])):
            loggingCritical("Error, the length of the different lists of qnbs are not valid in Fermions.",
                            verbose=RunSettings['vCritical'])
            exit()
        else:
            for part, qnbs in yamlSettings['Fermions'].items():
                tp1 = {}
                for igg, gg in enumerate(ListG):
                    tp1.update({gg: qnbs[igg + 1]})
                yamlSettings['Fermions'][part] = {'Gen': qnbs[0], 'Qnb': tp1}
    if 'RealScalars' in yamlSettings and yamlSettings['RealScalars'] != {} and type(
            yamlSettings['RealScalars'].values()[0]) == list:
        if not (all([len(pp) == len(yamlSettings['Groups']) for pp in yamlSettings['RealScalars'].values()])):
            loggingCritical("Error, the length of the different lists of qnbs are not valid in RealScalars.",
                            verbose=RunSettings['vCritical'])
            exit()
        else:
            for part, qnbs in yamlSettings['RealScalars'].items():
                tp1 = {}
                for igg, gg in enumerate(ListG):
                    tp1.update({gg: qnbs[igg]})
                yamlSettings['RealScalars'][part] = copy.deepcopy(tp1)
    if 'CplxScalars' in yamlSettings and yamlSettings['CplxScalars'] != {} and type(
            yamlSettings['CplxScalars'].values()[0]['Qnb']) == list:
        if not (all([len(pp['Qnb']) == len(yamlSettings['Groups']) for pp in yamlSettings['CplxScalars'].values()])):
            loggingCritical("Error, the length of the different lists of qnbs are not valid in CplxScalars.",
                            verbose=RunSettings['vCritical'])
            exit()
        else:
            for part, qnbs in yamlSettings['CplxScalars'].items():
                tp1 = {}
                for igg, gg in enumerate(ListG):
                    tp1.update({gg: qnbs['Qnb'][igg]})
                yamlSettings['CplxScalars'][part]['Qnb'] = copy.deepcopy(tp1)
    if 'Potential' in yamlSettings and yamlSettings['Potential'] != {}:
        Labels = {'QuarticTerms': 4, 'Yukawas': 3, 'TrilinearTerms': 3, 'ScalarMasses': 2, 'FermionMasses': 2}
        # There two things two do:1. deal with short hand notation 2. Deal with list of lists
        for lab, ilab in Labels.items():
            if lab in yamlSettings['Potential'] and yamlSettings['Potential'][lab] != {}:
                for ill1, ll1 in yamlSettings['Potential'][lab].items():
                    if type(ll1) == list and type(ll1[0]) == list:
                        # it is a list of list
                        if not (all([len(subll1) == ilab or len(subll1) == ilab + 1 for subll1 in ll1])):
                            loggingCritical(
                                "Error, the length of the list {} is incoherent should be {} or {} if norm is included".format(
                                    ll2, ilab, ilab + 1))
                            exit()
                        else:
                            nn = []
                            templ1 = []
                            for ll2 in ll1:
                                nn.append(1) if len(ll2) == ilab else nn.append(ll2[-1])
                                templ1.append(ll2[:ilab])
                            yamlSettings['Potential'][lab][ill1] = {'Fields': templ1, 'Norm': nn}
                    elif type(ll1) == list:
                        if not (len(ll1) == ilab or len(ll1) == ilab + 1):
                            loggingCritical(
                                "Error, the length of the list {} is incoherent should be {} or {} if norm is included".format(
                                    ll1, ilab, ilab + 1))
                            exit()
                        else:
                            nn = 1 if len(ll1) == ilab else ll1[-1]
                            yamlSettings['Potential'][lab][ill1] = {'Fields': ll1[:ilab], 'Norm': nn}
                    else:
                        # It is a regular dictionary we have to check that the norm is understood correctly i.e.
                        # 1) Only one number and Fields is a simple list ->Ok
                        # 2) Only one number and Fields is a list of lsit -> all the same norm
                        # 3) List of numbers and Fields is a list of list -> Ok
                        if 'Norm' in ll1:
                            if type(ll1['Fields'][0]) == list and (
                                            type(ll1['Norm']) == str or type(ll1['Norm']) == int):
                                # copy it this is case 2)
                                yamlSettings['Potential'][lab][ill1] = {'Fields': ll1['Fields'],
                                                                        'Norm': [ll1['Norm']] * len(ll1['Fields'])}
                            elif type(ll1['Fields'][0]) == list and type(ll1['Norm']) == list:
                                pass
                            elif type(ll1['Fields'][0]) == str and (
                                            type(ll1['Norm']) == str or type(ll1['Norm']) == int):
                                pass
                            else:
                                loggingCritical(
                                    "Error format for the Norm not understood. Should be either a list of same length as Fields or a single number applying to all of them: {}".format(
                                        ll1), verbose=True)
                                exit()
                        else:
                            loggingInfo("Warning no Norm specify for {}, setting it to `1`".format(ll1),
                                        verbose=RunSettings['vInfo'])
                            if type(ll1['Fields'][0]) == list:
                                yamlSettings['Potential'][lab][ill1] = {'Fields': ll1['Fields'],
                                                                        'Norm': [1] * len(ll1['Fields'])}
                            else:
                                yamlSettings['Potential'][lab][ill1] = {'Fields': ll1['Fields'], 'Norm': 1}
    # import the module
    loggingInfo("Importing modules ...", verbose=RunSettings['vInfo'])
    from ModelsClass import *
    from RGEsModule import *

    loggingInfo("\t\t...done", verbose=RunSettings['vInfo'])

    # Convert the quantum numbers into Rational of Sympy and the particles
    if 'Fermions' in yamlSettings:
        for key, value in yamlSettings['Fermions'].items():
            try:
                for GaugeGroup, qnb in value['Qnb'].items():
                    if type(qnb) == list:  # given in terms of the Dynkin label translate it into tuple
                        yamlSettings['Fermions'][key]['Qnb'][GaugeGroup] = tuple(qnb)
                    else:
                        try:
                            yamlSettings['Fermions'][key]['Qnb'][GaugeGroup] = Rational(qnb)
                        except ValueError as err:
                            yamlSettings['Fermions'][key]['Qnb'][GaugeGroup] = eval(qnb)
                    loggingInfo("Fermion {} with Qnb {} under {}".format(key, qnb, GaugeGroup))
            except ValueError as err:
                loggingCritical(
                    "Error during convertion of the quantum numbers : {} {} {}.".format(key, GaugeGroup, qnb),
                    verbose=RunSettings['vCritical'])
                exit()
            except AttributeError as err:
                loggingCritical("Error, the Qnb entry must be a dictionary", versbose=RunSettings['vCritical'])
                exit()
    if 'RealScalars' in yamlSettings:
        for key, value in yamlSettings['RealScalars'].items():
            try:
                for GaugeGroup, qnb in value.items():
                    if type(qnb) == list:  # given in terms of the Dynkin label translate it into tuple
                        yamlSettings['RealScalars'][key]['Qnb'][GaugeGroup] = tuple(qnb)
                    else:
                        yamlSettings['RealScalars'][key][GaugeGroup] = Rational(qnb)
                    loggingInfo("RealScalar {} with Qnb {} under {}".format(key, qnb, GaugeGroup))
            except ValueError as err:
                loggingCritical(
                    "Error during convertion of the quantum numbers : {} {} {}.".format(key, GaugeGroup, qnb),
                    verbose=RunSettings['vCritical'])
                exit()
            except AttributeError as err:
                loggingCritical("Error, the Qnb entry must be a dictionary", versbose=RunSettings['vCritical'])
                exit()
    if 'CplxScalars' in yamlSettings:
        for key, value in yamlSettings['CplxScalars'].items():
            try:
                for GaugeGroup, qnb in value['Qnb'].items():
                    if type(qnb) == list:  # given in terms of the Dynkin label translate it into tuple
                        yamlSettings['CplxScalars'][key]['Qnb'][GaugeGroup] = tuple(qnb)
                    else:
                        yamlSettings['CplxScalars'][key][GaugeGroup] = Rational(qnb)
                    loggingInfo("CplxScalar {} with Qnb {} under {}".format(key, qnb, GaugeGroup))
            except ValueError as err:
                loggingCritical(
                    "Error during convertion of the quantum numbers : {} {} {}.".format(key, GaugeGroup, qnb),
                    verbose=RunSettings['vCritical'])
                exit()
            except AttributeError as err:
                loggingCritical("Error, the Qnb entry must be a dictionary", versbose=RunSettings['vCritical'])
                exit()
        loggingInfo("Translating the quantum numbers into Rational of Sympy... done", verbose=RunSettings['vDebug'])

    # Kin mixing
    yamlSettings['KinMix'] = RunSettings['KinMix']
    # Create the instance Model
    model = Model(yamlSettings)
    if model.kinmixing:
        loggingInfo("Kinetic Mixing switched on, calculating kinetic mixing.",verbose=RunSettings['vInfo'])
    loggingCritical(
        "Creating the instance of the Model: {}, {}, by {}...done".format(yamlSettings['Name'], yamlSettings['Date'],
                                                                          yamlSettings['Author']),
        verbose=RunSettings['vCritical'])
    if RunSettings['SetGutNorm'] and len(model.UGaugeGroups) == 1:
        loggingInfo('Setting the Gut normalization for the U(1) group factor', verbose=RunSettings['vInfo'])
        RunSettings['SetGutNorm'] = True
    elif RunSettings['SetGutNorm'] and not (len(model.UGaugeGroups) == 1):
        loggingCritical('**Warning** the `SetGutNorm` switch ignored because multiple U(1) gauge groups found',
                        verbose=RunSettings['vCritical'])
        RunSettings['SetGutNorm'] = False

    # add the element in ToOnly to Only
    # IF only is defined then skip the other terms
    trans = {'QuarticTerms': model.LbdToCalculate,
             'ScalarsMasses': model.ScMToCalculate,
             'FermionMasses': model.FMToCalculate,
             'TrilinearTerms': model.TriToCalculate,
             'Yukawas': model.YukToCalculate,
             'ScalarAnomalous': model.ScalarAnomalousToCalculate,
             'FermionAnomalous': model.FermionAnomalousToCalculate
             }

    if RunSettings['Only'] != []:
        if type(RunSettings['Only']) == str:
            try:
                RunSettings['Only'] = eval(RunSettings['Only'])
            except:
                loggingInfo(
                    'WARNING: wrong format for the `Only` switch, take care of putting the `"` correctly, ignoring it',
                    verbose=True)
                RunSettings['Only'] = []
        todelete = []
        for valonly in RunSettings['Only']:
           for key,val in trans.items():
               if valonly in val.keys():
                   # we have identified which kind of coupling it is
                   for vv in val.keys():
                       if not vv in RunSettings['Only']:
                           trans[key].__delitem__(vv)
                           loggingDebug("\tIgnoring `{}` from `{}` according to `Only` settings".format(vv,key),verbose=True)


    # Check which terms are in the Lagrangian and if the corresponding rges have been calculated
    Translation = {'Yukawas': 'Yukawas', 'QuarticTerms': 'Quartic-Couplings', 'ScalarMasses': 'ScalarMass',
                   'FermionMasses': 'FermionMass', 'TrilinearTerms': 'Trilinear'}
    ListLagrangian = model.Potential.keys()
    StrucYuk = {}
    if RunSettings['Export']:
        Check = all(
            [True if RunSettings[Translation[el]] else False for el in ListLagrangian if model.Potential[el] != {}])
        Check = Check and RunSettings['Gauge-Couplings']
        if not (Check) and (not (RunSettings['All-Contributions'])) or RunSettings['Only'] != []:
            loggingCritical(
                "\n**WARNING** : options --Export cannot be carried out because all the required terms have not been set to True or Only is being used.\n Please use the --All-Contributions switch or set all the terms to True.\n",
                verbose=RunSettings['vCritical'])
            Check = False
        else:
            Check = True
    # check that there re no symbols in the definition of the generation
    CheckNumerics = True
    if 'Yukawas' in model.Potential:
        for key, val in model.Potential['Yukawas'].items():
            StrucYuk[key] = [model.Fermions[str(el.base)].Gen for el in val[0][0] if str(el.base) in model.Fermions]
            StrucYuk[key].append(
                [model.Scalars[str(el.base)].Gen for el in val[0][0] if str(el.base) in model.Scalars][0])
            # Check that there's only number here
            CheckNumerics = CheckNumerics and all([True if type(el) != Symbol else False for el in StrucYuk[key]])
    if 'FermionMasses' in model.Potential and CheckNumerics:
        for key, val in model.Potential['FermionMasses'].items():
            StrucYuk[key] = [model.Fermions[str(el.base)].Gen for el in val[0][0] if str(el.base) in model.Fermions]
            # Check that there's only number here
            CheckNumerics = CheckNumerics and all([True if type(el) != Symbol else False for el in StrucYuk[key]])
    if not (CheckNumerics) and RunSettings['Export']:
        loggingCritical(
            "\n**WARNING** : options --Export cannot be carried out because some of the generation indices are symbols.",
            verbose=RunSettings['vCritical'])
        Check = False

    # Start the actual calculation of the RGEs
    loggingInfo("Starting the Calculation...\n", verbose=RunSettings['vInfo'])
    ToCalculate = []
    if RunSettings['All-Contributions']:
        strListLagrangian = '\n\t\t\t\tGauge-Couplings'
        for el in ListLagrangian:
            strListLagrangian += '\n\t\t\t\t' + str(el)
        loggingInfo("\t\t... with option --All-Contributions including{} ".format(strListLagrangian),
                    verbose=RunSettings['vInfo'])
        ToCalculate = ['Gauge-Couplings'] + [Translation[el] for el in ListLagrangian] + ['ScalarAnomalous'] + ['FermionAnomalous']
    else:
        if RunSettings['Gauge-Couplings']:
            ToCalculate.append('Gauge-Couplings')
            loggingInfo("\t\t... with option --Gauge-Couplings", verbose=RunSettings['vInfo'])
        if RunSettings['Quartic-Couplings']:
            ToCalculate.append('Quartic-Couplings')
            loggingInfo("\t\t... with option --Quartic-Couplings", verbose=RunSettings['vInfo'])
        if RunSettings['Yukawas']:
            ToCalculate.append('Yukawas')
            loggingInfo("\t\t... with option --Yukawas", verbose=RunSettings['vInfo'])
        if RunSettings['ScalarMass']:
            ToCalculate.append('ScalarMass')
            loggingInfo("\t\t... with option --ScalarMass", verbose=RunSettings['vInfo'])
        if RunSettings['FermionMass']:
            ToCalculate.append('FermionMass')
            loggingInfo("\t\t... with option --FermionMass", verbose=RunSettings['vInfo'])
        if RunSettings['Trilinear']:
            ToCalculate.append('Trilinear')
            loggingInfo("\t\t... with option --Trilinear", verbose=RunSettings['vInfo'])
        if RunSettings['ScalarAnomalous']:
            if model.ScalarAnomalousToCalculate == {}:
                loggingCritical(
                    "\t\t**WARNING**, No Scalar Anomalous dimension found to calculate. Each anomalous is associated to a scalar mass term, make sure that the corresponding scalar mass terms are defined in the potential.",
                    verbose=True)
            else:
                ToCalculate.append('ScalarAnomalous')
                loggingInfo("\t\t... with option --ScalarAno", verbose=RunSettings['vInfo'])
        if RunSettings['FermionAnomalous']:
            if model.FermionAnomalousToCalculate == {}:
                loggingCritical(
                    "\t\t**WARNING**, No Fermion Anomalous dimension found to calculate. Each anomalous is associated to a fermion mass term, make sure that the corresponding fermion mass terms are defined in the potential.",
                    verbose=True)
                ToCalculate.append('FermionAnomalous')
            else:
                ToCalculate.append('FermionAnomalous')
                loggingInfo("\t\t... with option --FermionAno", verbose=RunSettings['vInfo'])


    # To store the Final Result
    FinalRGEdic = {}
    FinalRGE = []
    for rge in ToCalculate:
        temp = Translate(rge, model, RunSettings)
        FinalRGEdic[rge] = temp
        FinalRGE.append(temp)
    FinalRGEForTranslation = copy.deepcopy(FinalRGEdic)
    FinalRGEdic = FinalRGEdic.values()
    loggingCritical(" \t\t...Calculation done.", verbose=RunSettings['vCritical'])
    # Copy the settings in model for simplicity
    model.savesettingsforoutput.append(ToCalculate)
    model.savesettingsforoutput.append(StrucYuk)
    model.savesettingsforoutput.append(FinalRGEForTranslation)
    model.savesettingsforoutput.append(RunSettings)
    ###########################################################################################################################
    # Cleanup the RGE dictionary so that we can use the result more easily. This consists in defining a Function for each symbol
    # in classes
    ###########################################################################################################################
    FuncDefined = {}
    Subs = []
    for key, val in model.Classes.items():
        kkey = ''.join(reg.split('{|}', str(key)))
        kkeyf = Function(kkey, commutative=False)
        FuncDefined[key] = kkeyf, Symbol(kkey, commutative=False)
    # lets get the symbols to replace
    for yuk in model.ListYukawa + model.ListFM:
        # Get the rge
        rge = [el for el in FinalRGE if yuk in el]
        if rge == []:
            loggingDebug("Could not simplify the expression because all the necessary rges have not been calculated",
                         verbose=RunSettings['vDebug'])
        else:
            rge = rge[0]
            # Get the external symbols
            # two symbols
            extsymbs = rge[yuk].replace(model.Classes[yuk](p, q), Integer(1), map=True)[-1].keys()
            if extsymbs != []:
                extsymbs = extsymbs[0]
            else:
                # one symbols
                extsymbs = rge[yuk].replace(model.Classes[yuk](p), Integer(1), map=True)[-1].keys()
                if extsymbs != []:
                    extsymbs = extsymbs[0]
                    Subs.append((extsymbs, FuncDefined[yuk][0](*extsymbs.arg)))
            Subs.append((Symbol(yuk, commutative=False), FuncDefined[yuk][-1]))
    RGEs = []
    for rge in FinalRGE:
        tp = {}
        for key, val in rge.items():
            if key == 'abelian':
                tp[key] = [el.subs(Subs) for el in val]
            else:
                tp[key] = val.subs(Subs)
        RGEs.append(copy.deepcopy(tp))

    # Importing the requesting modules
    from ToMathematica import TranslateToMathematica
    from ExportBetaFunction import ExportBetaFunction, ExportBetaToCpp
    from latexRGEs import writeLatexOutput

    #############
    # LateX output
    #############

    # protect the export with try statement so that if it crashes it can continue
    
    try:
        if RunSettings['LatexOutput']:
            loggingInfo("\nWriting a latex file `{}` .".format(RunSettings['Results'] + '/' + RunSettings['LatexFile']),
                        verbose=RunSettings['vInfo'])
            # Call the writeLateX function
            writeLatexOutput(RunSettings, FinalRGE, model, ToCalculate)
    except:
        loggingCritical('"\nError in the LateX export, skipping. Contact the authors', verbose=True)
        pass

    ###################
    # Mathematica Output
    ###################
    try:
        if RunSettings['ToM']:
            loggingInfo(
                "Writing a txt to Mathematica file `{}`.".format(RunSettings['Results'] + '/' + RunSettings['ToMF']),
                verbose=RunSettings['vInfo'])
            AllMathematicaRGEs = []
            Dimension = []
            for iel, el in enumerate(FinalRGE):
                tp = TranslateToMathematica(el, '{}_{}'.format(RunSettings['Results'] + '/' + RunSettings['ToMF'],
                                                               ToCalculate[iel]), ToCalculate[iel], model)
                AllMathematicaRGEs.append(tp[0])
                Dimension.append(tp[1])
            f = open(RunSettings['Results'] + '/' + RunSettings['ToMF'].split('.txt')[0] + '_numerics.m', 'w')
            Dimension = 'DimensionParameter={{{}}};\n'.format(','.join(flatten(Dimension)))
            AllMathematicaRGEs = 'AllRGEs={{{}}}'.format(','.join(flatten(AllMathematicaRGEs)))
            f.write(Dimension)
            f.write('\n')
            f.write(AllMathematicaRGEs)
            f.close()
    except:
        loggingCritical('"\nError in the Mathematica export, skipping. Contact the authors', verbose=True)
        pass


    ##############
    # Pickle Output
    ##############
    try:
        if RunSettings['Pickle']:
            try:
                fPickle = open('{}'.format(RunSettings['Results'] + '/' + RunSettings['PickleFile']), 'w')
            except:
                loggingCritical('ERROR while opening the pickle output file `{}`.'.format(fPickle),
                                verbose=RunSettings['vCritical'])
            try:
                # create a string of the result
                if 'abelian' in RGEs[0]:
                    for iel, el in enumerate(flatten(model.UsectorMatrix)):
                        RGEs[0][str(el)] = RGEs[0]['abelian'][iel]
                    del (RGEs[0]['abelian'])

                strres = [str(el).replace('Dagger', 'adjoint') for el in RGEs]
                if model.kinmixing:
                    gaugegroups = [str(el) for el in model.UsectorMatrix] + [str(el[1].g) for el in model.NonUGaugeGroups]
                else:
                    gaugegroups = [str(el[-1]) for el in model.ListGroups]
                pudb.set_trace()
                strsettings = ([str(el) for el in model.ListYukawa],
                               [str(el) for el in model.ListFM],
                               [str(el) for el in model.ListTri],
                               [str(el) for el in model.ListLbd],
                               [str(el) for el in model.ListScM],
                               [el for el in gaugegroups],
                               [str(el[0]) for el in model.GaugeGroups],
                               [str(el) for el in model.Particles.keys()],
                               [str(el) + '_f' for el in model.Particles.keys()],
                               [str(el) + '1_f' for el in model.Particles.keys()],  # for the indices of the generations
                               [str(el) for el in model.ScalarAnomalousToCalculate.keys()],
                               [str(el) for el in model.FermionAnomalousToCalculate.keys()],
                               ['xi', 'g_SUndum']
                               )
                # Fix for the Sudummy gauge group
                if 'SUndum' in strsettings[-4]:
                    strsettings[-5].append('g_SUndum')
                pickle.dump([strres, strsettings, model.saveSettings, 'Generated by PyR@TE for {}.'.format(model._Name)], fPickle)
                fPickle.close()
                loggingInfo("\nWriting a pickle output file `{}`.".format(RunSettings['PickleFile']),
                            verbose=RunSettings['vInfo'])
            except:
                loggingCritical('ERROR while storing the pickle output file.', verbose=RunSettings['vCritical'])
    except:
        loggingCritical('"\nError in the Pickle export, skipping. Contact the authors', verbose=True)
        pass

    ########################
    # Python numerical Output
    ########################
    try:
        if RunSettings['Export'] and Check:
            loggingInfo("Exporting the results to a python function...", verbose=RunSettings['vInfo'])
            pythonoutput = ExportBetaFunction(model, FinalRGEForTranslation, RunSettings, StrucYuk)
            loggingInfo("\t\t...done", verbose=RunSettings['vInfo'])
            loggingInfo("Exporting the results to a c++ function...", verbose=RunSettings['vInfo'])
            ExportBetaToCpp(pythonoutput, 'BetaFunction{}.cpp'.format(model._Name))
            loggingInfo("\t\t...done", verbose=RunSettings['vInfo'])
    except:
        loggingCritical('"\nError in the numerical export, skipping. Contact the authors', verbose=True)
        pass
    model.idb._update_db()
    loggingCritical('End of the run.', verbose=RunSettings['vCritical'])
    # mv the logging into results
    os.system('mv {} {}'.format(RunSettings['LogFile'], RunSettings['Results']))

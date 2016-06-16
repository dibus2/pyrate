"""

This is a set of wrapping function for the ipython notebook.
The first one is the LateX output then we have the pickle and Mathematica and finally the Load function 

"""
try:
    import sys
    import os

    wd = os.getcwd()
    sys.path.append(wd + '/src/Output')
    sys.path.append(wd + '/src/GroupTheory/')
    sys.path.append(wd + '/src/Core')
    import string
    import time as tm
    import re as r
except ImportError:
    exit("error while loading the modules in Toolbox")
from Logging import *
from latexRGEs import *
from ToMathematica import TranslateToMathematica, ToMathematicaNotation, findclosingbracket
from loadingEqs import Load
from sympy import transpose
from RGEsmathModule import trace, SP, MatM, Tr, adjoint


###########################################
# LateX Output 
###########################################

def ExportToLateX(Name, Xpr, model):
    # open the file
    try:
        fLatex = open(Name, 'w')
    except:
        print("error, impossible to open the file {}".format(Name))
    print"Exporting to LateX..."
    LatexText = r"""
		\documentclass[12pt,letter,fleqn]{article}
		\pagestyle{empty}
		\usepackage[latin1]{inputenc}
		\usepackage[dvips,landscape,top=1cm,nohead,nofoot]{geometry}
		\usepackage{amsmath}
		\usepackage{amsfonts}
		\usepackage{amssymb}
		\usepackage{scrtime}
		\setlength{\parindent}{0pt}
		\newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}
		\newcommand{\conj}[1]{\ensuremath{#1^*}}
		\begin{document}
		\begin{flushright}\today \quad \thistime\end{flushright}"""
    LatexText += r"""\begin{{center}} model: {} \end{{center}}""".format(model._Name)
    counter = 0
    Translation = {
        'Gauge-Couplings': "writeGaugeCouplings(LatexText,Xpr,model,model.savesettingsforoutput[-1])",
        'Quartic-Couplings': "writeQuarticCouplings(LatexText,Xpr,model,model.savesettingsforoutput[-1],counter)",
        'Yukawas': "writeYukawas(LatexText,Xpr,model,model.savesettingsforoutput[-1],counter)",
        'ScalarMass': "writeScalarMass(LatexText,Xpr,model,model.savesettingsforoutput[-1],counter)",
        'FermionMass': "writeFM(LatexText,Xpr,model,model.savesettingsforoutput[-1],counter)",
        'Trilinear': "writeTrilinear(LatexText,Xpr,model,model.savesettingsforoutput[-1],counter)"
    }
    for el in model.savesettingsforoutput[0]:
        LatexText = eval(Translation[el])
        counter += 1
    LatexText += '\end{document}'
    fLatex.write(LatexText)
    fLatex.close()
    print"\t\t...done"


def ExportToMathematica(Name, Xpr, model):
    # Loop other the results
    print "writing a Mathematica txt file..."
    AllMathematicaRGEs = []
    Dimension = []
    for iel, el in enumerate(Xpr):
        tp = TranslateToMathematica(el, '{}_{}'.format(Name, model.savesettingsforoutput[0][iel]),
                                    model.savesettingsforoutput[0][iel], model)
        AllMathematicaRGEs.append(tp[0])
        Dimension.append(tp[1])
    f = open(Name + '_numerics.m', 'w')
    Dimension = 'DimensionParameter={{{}}};\n'.format(','.join(flatten(Dimension)))
    AllMathematicaRGEs = 'AllRGEs={{{}}}'.format(','.join(flatten(AllMathematicaRGEs)))
    f.write(Dimension)
    f.write('\n')
    f.write(AllMathematicaRGEs)
    f.close()
    print"\t\t...done"


def ExportToPickle(Name, Xpr, model, description=''):
    fPickle = open(Name, 'w')
    # create a string of the result
    print "Dumping the results in a pickle file..."
    if 'abelian' in Xpr:
        for iel, el in enumerate(flatten(model.UsectorMatrix)):
            RGEs[0][str(el)] = RGEs[0]['abelian'][iel]
        del (RGEs[0]['abelian'])
    if type(Xpr) == dict:
        strres = [str(el).replace('Dagger', 'adjoint') for el in [Xpr]]
    else:
        strres = [str(el).replace('Dagger', 'adjoint') for el in Xpr]
    if model.kinmixing:
        gaugegroups = [str(el) for el in model.UsectorMatrix] + [str(el[1].g) for el in model.NonUGaugeGroups]
    else:
        gaugegroups = [str(el[-1]) for el in model.ListGroups]
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
                   ['xi', 'g_SUndum']
                   )
    # Fix for the Sudummy gauge group
    if 'SUndum' in strsettings[-4]:
        strsettings[-5].append('g_SUndum')

    pickle.dump([strres, strsettings, model.saveSettings, description], fPickle)
    fPickle.close()
    print"\t\t...done"


def loadmodel(name):
    # open the file
    Of = open('{}'.format(name), 'r')
    strres, strsettings, modelsettings, description = pickle.load(Of)
    # we first declare all the variables and we also have to declare the classes to make sense of the different expressions
    Replace = []
    for ill, ll in enumerate(strsettings):
        # simplify the symbol to get rid of the LateX structure
        lln1 = [(iel, ''.join(reg.split('{(.*)}', el))) for iel, el in enumerate(ll)]
        lln1 = [(iel, ''.join(reg.split(r'\\', el))) for iel, el in enumerate([elem[1] for elem in lln1])]
        lln = lln1
        Replace.append([(ll[iel], el) for iel, el in lln])
        # Declare all the symbols
        ll = ll + [el[1] for el in lln]
        for symb in ll:
            if ill != 0 and ill != 1:
                globals()[symb] = Symbol(symb)
            else:
                globals()[symb] = Symbol(symb, commutative=False)
    Replace = sum(Replace, [])
    for rr in Replace:
        strres = [el.replace(*rr) for el in strres]
    res = [eval(el) for el in strres]
    # translate the keys into the symbols
    Nres = {}
    pudb.set_trace()
    for elem in res:
        for key, vals in elem.items():
            Nres[key] = vals
    return Nres, description, modelsettings


def getoneloop(llist):
    if type(llist) == dict:
        return [(el * (4 * pi) ** 2).expand().subs(1 / pi ** 2, 0) for key, el in llist.items()]
    elif type(llist) == list:
        return [(el * (4 * pi) ** 2).expand().subs(1 / pi ** 2, 0) for el in llist]
    elif type(llist) == Mul or type(llist) == Add:
        return (llist * (4 * pi) ** 2).expand().subs(1 / pi ** 2, 0)
    else:
        print("wrong format use , dict, list, Mul or Add")
        return llist


def gettwoloop(llist):
    if type(llist) == dict:
        return [(el * (4 * pi) ** 4).expand().subs(pi ** 2, 0) for key, el in llist.items()]
    elif type(llist) == list:
        return [(el * (4 * pi) ** 4).expand().subs(pi ** 2, 0) for el in llist]
    elif type(llist) == Mul or type(llist) == Add:
        return (llist * (4 * pi) ** 4).expand().subs(pi ** 2, 0)
    else:
        print("wrong format use , dict, list, Mul or Add")
        return llist


def getdiff(a, b, subs=[]):
    assert len(a) == len(b)
    if subs != []:
        return [(el - b[i]).subs(*subs) if type(el - b[i]) != int else el - b[i] for i, el in enumerate(a)]
    else:
        tt = [el - b[i] for i, el in enumerate(a)]
        return tt


def CompareModels(model1, model2, subrules=[], display=False):
    if type(model1) == str and type(model2) == str:
        # load the models
        model1, description1 = loadmodel(model1)
        model2, description2 = loadmodel(model2)
    Diff = []
    if display:
        print "Difference of model1 - model2:"
    for key, rges in model2.items():
        # 1Loop
        if type(rges) != dict:
            rges = {'dummy': rges}
        model21L = getoneloop(rges)
        model22L = gettwoloop(rges)
        if key in model1:
            if type(model1[key]) != dict:
                model1[key] = {'dummy': model1[key]}
            model11L = getoneloop(model1[key])
            model12L = gettwoloop(model1[key])
        else:
            model11L = model12L = [0]
        diff1L = getdiff(model11L, model21L, subs=subrules)
        diff2L = getdiff(model12L, model22L, subs=subrules)
        Diff.append([key, diff1L, diff2L])
        if display:
            print "diff 1L {} : {}".format(key, diff1L)
        if model12L != [0] and model22L != [0]:
            if display:
                print "diff 2L {} : {}".format(key, diff2L)
    return Diff


###########################
# Export to a numerical file
###########################

def ExportBetaFunction(name, model):
    print"Exporting to python function..."
    cpStrucYuk = copy.deepcopy(model.savesettingsforoutput[1])
    settings = copy.deepcopy(model.savesettingsforoutput[-1])
    FinalExpr = copy.deepcopy(model.savesettingsforoutput[2])
    # kappa
    kappa = (4 * pi) ** 2
    # WARNING the scalars arre summed to all components of teh matrices
    date = tm.localtime()
    date = "{}-{}-{}\t {}:{}:{}".format(date[0], date[1], date[2], date[3], date[4], date[5])
    CodeStr = '\"\"\"\nThis Code is generated automatically by pyR@TE, {} for the model : {}\n\"\"\"\n'.format(
        date, model._Name)
    CodeStr += "from __future__ import division\nimport numpy as np\ndef beta_function_{}(t,y,Assumptions={{'two-loop':False,'diag': True}}):\n\t#The Assumption can be used to incorporate some switch in the calculation\n".format(
        model._Name)
    # We need an identity for each yukawa in principle because they could have different shape
    CodeStr += "\tkappa = 1./(16*np.pi**2)\n"
    # Simplify the keys of cpStrucYuk
    maps = {}
    Ids = {}
    for ikk, (key, val) in enumerate(cpStrucYuk.items()):
        maps[key] = ''.join(reg.split('\\\(.*)', ''.join(reg.split('{|}|_', key))))
        Ids[maps[key]] = "Id{}".format(ikk)
        if max(val) != 1:
            CodeStr += "\tId{} = np.eye({})\n".format(ikk, max(val))
        else:
            CodeStr += "\tId{} = 1\n".format(ikk)
    # We've got to get all the different RGEs and create variables
    NbRGEs = 0
    # List of all the symbols
    ListSymbs = [el[1].g for el in model.GaugeGroups] + model.ListLbd + model.ListScM + model.ListTri
    for iel, el in enumerate(ListSymbs):
        if len(reg.split('{(.*)}', str(el))) == 3 or len(reg.split('\\\(.*)', str(el))) == 3:
            ListSymbs[iel] = ''.join(reg.split('\\\(.*)', ''.join(reg.split('{|}|_', str(el)))))
            maps[el] = ListSymbs[iel]
        else:
            ListSymbs[iel] = str(el).replace('_', '').replace('\\', '').replace('\\\\', '').replace('^','').replace(' ', '')
            maps[el] = ListSymbs[iel]
    # First the scalar quatities
    NbRGEs += len(ListSymbs)
    # Now the non scalars quantities go through all of them and count the multiplicity due to the number of Generations
    if 'Yukawas' in model.Potential:
        ToLoop = {}
        ToLoop.update(model.Potential['Yukawas'])
        if 'FermionMasses' in model.Potential:
            ToLoop.update(model.Potential['FermionMasses'])
        for key, val in ToLoop.items():
            # Is it a squared matrix ? without the scalars
            if len(cpStrucYuk[key]) == 3:
                if cpStrucYuk[key][-1] != 1:
                    exit("Index families on higgs are not taken into account")
                else:
                    cpStrucYuk[key] = cpStrucYuk[key][:-1]
            # it can be scalar
            if cpStrucYuk[key] == [1, 1]:
                NbYuk = 1
                cpStrucYuk[key].append('scalar')
            else:
                NbYuk = cpStrucYuk[key][1] * cpStrucYuk[key][0]
                cpStrucYuk[key].append('matrix')
            NbRGEs += NbYuk
    # Count the rges that we have written
    iwritten = 0
    # declare the matrices
    Mapping = {}

    ##########################################################################################
    # Aux function
    def declarematrices(llist, iwritten, Mapping, cpStrucYuk, CodeStr):
        for ii, val in llist:
            if cpStrucYuk[ii][-1] == 'matrix':
                CodeStr += "\t{} = np.matrix(np.zeros(({},{})))\n".format(maps[ii], cpStrucYuk[ii][0],
                                                                          cpStrucYuk[ii][1])
                Mapping[ii] = (iwritten, iwritten + cpStrucYuk[ii][0] * cpStrucYuk[ii][1] - 1)
                # fill up the newly constructed matrix
                for j in range(cpStrucYuk[ii][0]):
                    for k in range(cpStrucYuk[ii][1]):
                        CodeStr += "\t{0}[{1},{2}] = y[{3}]\n".format(maps[ii], j, k, iwritten)
                        iwritten += 1
            elif cpStrucYuk[ii][-1] == 'scalar':
                Mapping[ii] = (iwritten, iwritten)
                # fill up the newly constructed matrix
                CodeStr += "\t{} = y[{}]\n".format(maps[ii], iwritten)
                iwritten += 1
        return Mapping, iwritten, CodeStr

    ##########################################################################################
    if 'Yukawas' in FinalExpr.keys():
        Mapping, iwritten, CodeStr = declarematrices(FinalExpr['Yukawas'].items(), iwritten, Mapping, cpStrucYuk,
                                                     CodeStr)
    if 'FermionMass' in FinalExpr.keys():
        Mapping, iwritten, CodeStr = declarematrices(FinalExpr['FermionMass'].items(), iwritten, Mapping, cpStrucYuk,
                                                     CodeStr)
    # Do the mapping once and for all
    for el in ListSymbs:
        Mapping[el] = iwritten
        iwritten += 1
    CodeStr += "\t###############\n\t#beta functions\n\t###############\n"
    CodeStr += "\t#This is the mapping between rges and y[i]\n\t#{}\n".format(Mapping)
    temp1L = temp2L = ''
    for ii, val in FinalExpr.items():
        if ii != 'ScalarAnomalous' and ii != 'FermionAnomalous':
            for label, expr in val.items():
                # deal with kinetic mixing
                if label == 'abelian':
                    # need to iterate over the size of the matrix
                    row, col = model.UsectorMatrix.shape
                    for irow in range(row * col):
                        label = str(model.UsectorMatrix[irow])
                        temp1L +='\tb{} = ({})*kappa\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa * expr[irow]).expand().subs(1 / pi ** 2, 0), ListSymbs, maps[label], Ids,
                                                model, Mapping, isScalar=True)
                        )
                        temp2L += '\t\tb{0} = b{0} + ({1})*kappa**2\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa ** 2 * expr[irow]).expand().subs(pi, 0), ListSymbs, maps[label], Ids, model,
                                                Mapping, isScalar=True)
                        )
                else:
                    if ii == 'Gauge-Couplings':
                        label = model.GetGroupFromName[label].g
                    # Get the 1loop contribution
                    if ii != 'Yukawas' and ii != 'FermionMass':
                        temp1L += '\tb{} = ({})*kappa\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa * expr).expand().subs(1 / pi ** 2, 0), ListSymbs, maps[label], Ids,
                                                model, Mapping, isScalar=True))
                        # Get the 2loop contribution
                        temp2L += '\t\tb{0} = b{0} + ({1})*kappa**2\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa ** 2 * expr).expand().subs(pi, 0), ListSymbs, maps[label], Ids, model,
                                                Mapping, isScalar=True))
                    else:
                        temp1L += '\tbeta{} = ({})*kappa\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa * expr).expand().subs(1 / pi ** 2, 0), ListSymbs, maps[label], Ids,
                                                model, Mapping, isScalar=False))
                        # Get the 2loop contribution
                        temp2L += '\t\tbeta{0} = beta{0} + ({1})*kappa**2\n'.format(
                            maps[label],
                            TranslateToNumerics((kappa ** 2 * expr).expand().subs(pi, 0), ListSymbs, maps[label], Ids, model,
                                                Mapping, isScalar=False))
    # We regroup all the 1 loop and 2 loop contributions to avoid multiple if statements
    CodeStr += temp1L
    CodeStr += "\tif Assumptions['two-loop']:\n{}".format(temp2L)

    # Create the return and write the file
    def Sort(x):
        if type(x[1]) == int:
            return x[1]
        else:
            return x[1][0]

    SortMapping = sorted(Mapping.items(), key=Sort)
    Return, ReturnFull = "", ""
    Tags, TagsFull = [], []
    for el in SortMapping:
        if type(el[1]) != int and cpStrucYuk[el[0]][-1] == 'matrix':
            for ii in range(cpStrucYuk[el[0]][0]):
                for jj in range(cpStrucYuk[el[0]][1]):
                    if ii == jj:
                        Return += 'beta{}[{},{}], '.format(maps[el[0]], ii, ii)
                        Tags.append('\"beta{}{}\"'.format(maps[el[0]], ii, ii))
                    else:
                        # If it is a vector all component should be returned
                        if cpStrucYuk[el[0]][1] == 1:  # vector
                            Return += 'beta{}[{},{}], '.format(maps[el[0]], ii, jj)
                            Tags.append('\"beta{}[{},{}]\"'.format(maps[el[0]], ii, jj))
                        else:
                            Return += '0,'
                            Tags.append('0')
                    ReturnFull += 'beta{}[{},{}], '.format(maps[el[0]], ii, jj)
                    TagsFull.append('\"beta{}[{},{}]\"'.format(maps[el[0]], ii, jj))
        elif type(el[1]) != int and cpStrucYuk[el[0]][-1] == 'scalar':
            Return += 'beta{}, '.format(maps[el[0]])
            ReturnFull += 'beta{}, '.format(maps[el[0]])
            TagsFull.append('\"beta{}\"'.format(maps[el[0]]))
        elif type(el[1]) == int:
            Return += 'b{}, '.format(el[0])
            ReturnFull += 'b{}, '.format(el[0])
            TagsFull.append('\"b{}\"'.format(el[0]))
        else:
            exit("Error this is not implemented, contact the authors")
    CodeStr += "\tif Assumptions['diag'] :\n\t\treturn [{}]\n\telse :\n\t\treturn [{}]".format(Return[:-2], ReturnFull[
                                                                                                            :-2])  # remove the last coma
    # transform Tag
    Tags = eval('[{}]'.format(','.join(Tags)))
    TagsFull = eval('[{}]'.format(','.join(TagsFull)))
    CodeStr += "\n\n\n\n#rge = RGE(beta_function_{},{},labels={})\n".format(model._Name, NbRGEs, TagsFull)
    # Check that there are no points and that it ends with a .py
    spltName = name.split('.')
    if len(spltName) > 2:
        Name = spltName[0]
        for ii in range(1, len(spltName) - 1):
            Name = Name + spltName[ii]
        if spltName[-1] == 'py':
            Name = Name + '.py'
        else:  # we add the rest and add the dot py by force
            Name = Name + spltName[-1] + '.py'
    elif len(spltName) == 1:
        Name = spltName[0] + '.py'
    elif len(spltName) == 2 and spltName[-1] != 'py':
        Name = spltName[0] + spltName[1] + '.py'
    else:
        Name = spltName[0] + '.' + spltName[1]
    f = open(name, 'w')
    f.write(CodeStr)
    f.close()

    #################################
    # End of the beta_function.py file
    #################################
    # Now let's create an example file in which we demonstrate how to use the above file in combination with the RGEClass.py
    # Header
    CodeStrI = '\"\"\"\nThis Code is generated automatically by pyR@TE, {} for the model : {}\nIt demonstrates how to use the beta_function file produced by pyR@TE using the -e option to solve the rges using the RGEsclass.py module\"\"\"\n'.format(
        date, model._Name)
    CodeStrI += "#!/usr/bin/env python\nimport sys\n#This line should be modified if you move the result file\nsys.path.append('{}')\n#Same as above.The RGEsclass.py should be accessible\nsys.path.append('{}')\ntry :\n\tfrom RGEclass import *\nexcept ImportError :\n\texit('RGEclass.py not found.')\ntry :\n\timport scipy.interpolate as scp\nexcept ImportError :\n\texit('error while loading scipy')\nfrom {} import *\n".format(
        settings['Results'], '{}/src/Output'.format(wd), Name.split('.')[0])
    CodeStrI += "#In Order to use your result you have to create an instance of the RGE class following this\nrge = RGE(beta_function_{},{},labels={})\n#Note here that the labels are just used to label the results but should match the order in which the RGEs are solved.\n".format(
        model._Name, NbRGEs, TagsFull)
    CodeStrI += "\"\"\"\nOne can access the physical parameters as the top mass, the up quark yukawas, Z mass ...\nSince they are all attributes of the class object we have just created. The list is the following :\n\tPi = np.pi\n\talpha = 1./128.91\n\talphaS = 0.1184\n\tsinthetasq = 0.2316 #Sin squared\n\tMz = 91.1876 #in GeV\n\tMw = 80.399 # in GeV\n\tGf = 1.16637e-5 #in GeV ^-2\n\tCF = 4./3. #Casimir of the SU3 color\n\tNc = 3.\n\tQtop = 2./3.#top quark chargs\n\tvev = 246.22 # From (sqrt(2)Gf)^-1/2\n\tMt = 173.3 # the pole mass of the top\n\tU = [2.4e-3,1.27,165.3308]#mu,mc,mt (MSbar)mt=(Mt/(1+4alphas/3Pi))\n\tD = [4.75e-3,104e-3,4.19]#md,ms,mb\n\tL = [0.511e-3,105.66e-3,1.777]\n\tyU = [sqrt(2)/self.vev * el for el in self.U]\n\tyD = [sqrt(2)/self.vev * el for el in self.D]\n\tyL = [sqrt(2)/self.vev * el for el in self.L] \n\tg10 = sqrt(4*self.Pi * self.alpha)/sqrt(1-self.sinthetasq)#Initial value\n\tg20 = sqrt(4*self.Pi * self.alpha)/sqrt(self.sinthetasq)#Initial value\n\tg30 =  sqrt(4*self.Pi * self.alphaS)#Initial value\n\n\nTo modify one of them one would do :\n\"\"\"\nrge.Mt = 171.0\n"
    CodeStrI += "#Now one should specify the initial values for the {0} parameters of the rges\n########################\n#USER INPUT REQUIRED HERE\n########################\nY0 = [0.0]*{0}\n#assign this list to the initial values of the rges\nrge.Y0 = Y0".format(
        NbRGEs)
    CodeStrI += "\n\"\"\"the translation is the following : {}\"\"\"\n".format(Mapping)
    CodeStrI += "#Define the t values, using numpy: tmin,tmax,tstep.t=log_10(Scale in GeV)\ntmin,tmax,tstep = np.log(rge.Mz)/np.log(10),19,0.1\n"
    CodeStrI += "#Let's define the assumptions that can be passed through this dictionnary. If you want to add other than just two-loop you have to do the corresponding modification in {}.py\nAssumptions={{'two-loop':False,'diag': True}}\n".format(
        settings['ExportFile'])
    CodeStrI += "#Now we can finally solve the rges for the required values\nrge.solve_rges(tmin,tmax,tstep,Assumptions)\n"
    CodeStrI += "#Finally if one wants to access the result it can use the labels define above for an example see below. Note that the 't' entry is defined by default.\nfrom pylab import *\n#Some settings\nrc('axes',linewidth=2)\ngrid(True)\nfigprops = dict(figsize=(6., 5.))\nfigure(1,**figprops)\nadjustprops = dict(top=0.95, bottom=0.15, left=0.18, right=0.95)\nsubplots_adjust(**adjustprops)\nplot(rge.Sol['t'],rge.Sol['{3}'],label='{3} of {0} model')\nplot(rge.Sol['t'],rge.Sol['{1}'],label='{1} of {0} model')\nplot(rge.Sol['t'],rge.Sol['{2}'],label='{2} of {0} model')\n#Some more settings\nxlabel(r'$t$',fontsize =18)\nylabel('{1},{2} and {3} params of {0}')\nleg=gca().legend(ncol=3,loc=3,shadow=True,numpoints=1)\nshow()\n".format(
        model._Name, TagsFull[0], TagsFull[1], TagsFull[2])
    f = open(name.split('.py')[0] + '_' + 'SolveRGEs.py', 'w')
    f.write(CodeStrI)
    f.close()
    print("\t\t...done")


def ToNumpy(func, xpr, Mat=False):
    pos0 = xpr.find('{}('.format(func)) + len(func)
    lstring = list(xpr)
    # and now find the matching round bracket
    pos1 = findclosingbracket(lstring, pos0 + 1)
    # change all the comas in between the two brackets
    while True:
        tpfind = xpr[pos0:].find(',')
        pcoma = pos0 + tpfind if tpfind != -1 else len(lstring) - 1
        if pcoma < pos1:
            lstring[pcoma] = '*'
            pos0 = pcoma + 1
        else:
            break
    if not (Mat):
        return ''.join(lstring).replace(func, 'np.{}'.format(func[1:]), 1)
    else:
        if func == 'SScalarProd':
            xpr = ''.join(lstring)
            src = reg.search('SScalarProd\(([^\*]*)\*([^\)]*)\)', xpr)
            if src != None:
                xpr = xpr.replace('SScalarProd({}*{})'.format(src.group(1), src.group(2)),
                                  'np.trace(np.transpose({})*{})'.format(src.group(1), src.group(2)))
            else:
                exit("error while translating the to numpy notation")
            return xpr
        else:
            return ''.join(lstring).replace(func, '', 1)


def TAdj(xpr):
    pos0 = xpr.find('Adj(')
    # replace only the first occurence
    xpr = xpr.replace('Adj', 'np.transpose(np.conjugate', 1)
    lstring = list(xpr)
    pos1 = findclosingbracket(lstring, pos0 + len('np.transpose(np.conjugate('))
    # add a bracket to close the additional one we introduced above
    lstring.insert(pos1 + 1, ')')
    return ''.join(lstring)


def RemoveIndices(xpr, el, Yuk=True):
    # The flag is here to deal with the MatMul
    if Yuk:
        # Removes the indices of a Yuk Symbol
        while el + '[' in xpr:
            pos0 = xpr.find(el + '[') + len(el)
            lstring = list(xpr)
            # remove the squar bracket
            lstring[pos0] = ''
            pos1 = xpr[pos0 + 1:].find(']') + pos0 + 1
            # remove indices
            lstring = lstring[:pos0] + lstring[pos1 + 1:]
            xpr = ''.join(lstring)
    else:
        while 'MatMul(' in xpr:
            pos0 = xpr.find('MatMul(') + len('MatMul(')
            lstring = list(xpr)
            lstring[pos0] = ''
            pos1 = findclosingbracket(lstring, pos0 + 1)
            if lstring[pos1 + 1] == ',':
                pos2 = xpr[pos1 + 1:].find(')') + pos1 + 1
                lstring = lstring[:pos1 + 1] + lstring[pos2 + 1:]
            elif lstring[pos1 + 1] == ')' and lstring[pos1] == ')':
                lstring[pos1] = ''
            else:
                loggingDebug("Warning case not taken into account in RemoveIndices")
                pass
            xpr = ''.join(lstring)
            xpr = xpr.replace('MatMul', 'matMul', 1)
    return xpr


def TranslateToNumerics(expression, ListSymbs, label, Ids, model, Mapping, Matrices=True, isScalar=False):
    """Implements all the modification to the expression to be evaluated numerically and return a string"""
    # Use the ToMathematicaNotation to clean up the result a bit
    StrXpr = ToMathematicaNotation(expression, model, FlagSquare=False)  # to avoid the '('->'['
    if not (isScalar):
        StrXpr = StrXpr.replace('trace', '{}*ttrace'.format(Ids[label])).replace('conj', 'np.conjugate').replace('Tp',
                                                                                                                 'np.transpose').replace(
            'ScalarProd', '{}*SScalarProd'.format(Ids[label]))
    else:
        StrXpr = StrXpr.replace('trace', 'ttrace').replace('conj', 'np.conjugate').replace('Tp',
                                                                                           'np.transpose').replace(
            'ScalarProd', 'SScalarProd')
    while 'ttrace(' in StrXpr:
        StrXpr = ToNumpy('ttrace', StrXpr)
    while 'Adj(' in StrXpr:
        StrXpr = TAdj(StrXpr)
    # Gauge-Couplings,Quartic,Trilinear,scalar masses
    for el in ListSymbs:
        StrXpr = StrXpr.replace(str(el), 'y[{}]'.format(Mapping[el]))
    # Now let's deal with the Yukawas and the indices
    if Matrices:  # result in matrix form
        # Remove all the indices
        for el in model.ListYukawa:
            el = ''.join(reg.split('\\\(.*)', ''.join(reg.split('{|}|_', el))))
            StrXpr = RemoveIndices(StrXpr, el, Yuk=True)
        StrXpr = RemoveIndices(StrXpr, 'MatMul', Yuk=False)
        while 'matMul(' in StrXpr:
            StrXpr = ToNumpy('matMul', StrXpr,
                             Mat=True)  # The Mat is used to remove the tag of the function like ScalarProd -> ''
        while 'SScalarProd(' in StrXpr:
            StrXpr = ToNumpy('SScalarProd', StrXpr, Mat=True)
    return StrXpr


#########################################
# This part is not part of the public code
#########################################

# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def filtertraceMatM(xpr, dic):
    """
    check if there is a zero matrix in dic and return zero if it is the case.
    """
    if any([val == 0 for key, val in dic.items()]):
        return Integer(0)


def SimplifytraceMatM(xpr):
    """
    remove all the traces and MatM, SP that contain a zero matrix in.
    """
	# terms with two or four objects
    p, q, r, s, t, u, v, w = map(Wild, ['p', 'q', 'r', 's', 't', 'u', 'v', 'w'])
    xpr = xpr.replace(lambda xpr: xpr.match(trace(p, q)) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(trace(p, q))))
    xpr = xpr.replace(lambda xpr: xpr.match(trace(p, q, r, s)) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(trace(p, q, r, s))))
    xpr = xpr.replace(lambda xpr: xpr.match(MatM((p, q))) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(MatM((p, q)))))
    xpr = xpr.replace(lambda xpr: xpr.match(MatM((p, q, r, s))) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(MatM((p, q, r, s)))))
    xpr = xpr.replace(lambda xpr: xpr.match(MatM((p, q), r, s)) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(MatM((p, q), r, s))))
    xpr = xpr.replace(lambda xpr: xpr.match(MatM((p, q, u, v), r, s)) != None,
                      lambda xpr: filtertraceMatM(xpr, xpr.match(MatM((p, q, u, v), r, s))))
    xpr = xpr.replace(lambda xpr: xpr.match(SP(p, q)) != None, lambda xpr: filtertraceMatM(xpr, xpr.match(SP(p, q))))
    return xpr


def settozero(xpr, symbs, com=False):
    """
    put symbols to zeros.
    """
    if not (com):
        DeclaredSymbols = [Symbol('{}'.format(el), commutative=False) for el in symbs]
    else:
        DeclaredSymbols = [Symbol('{}'.format(el), commutative=True) for el in symbs]
    Subs = tuple([(el, 0) for el in DeclaredSymbols])
    return SimplifytraceMatM(xpr.subs(Subs))

# <codecell>

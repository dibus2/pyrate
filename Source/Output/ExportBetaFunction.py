try :
        import sys
        import os
        import platform
        import string as st
        cwd = os.getcwd()
        sys.path.append('../GroupTheory')
        from sympy import latex
except ImportError :
        exit('error while loading modules in latexRGEs.py')
from Contraction import *
from ToMathematica import ToMathematicaNotation,findclosingbracket

def ExportBetaFunction(model,FinalExpr,settings,StrucYuk):
        #kappa 
        kappa = (4*pi)**2
        #WARNING the scalars arre summed to all components of teh matrices
        date = time.localtime()
        date ="{}-{}-{}\t {}:{}:{}".format(date[0],date[1],date[2],date[3],date[4],date[5])
        CodeStr = '\"\"\"\nThis Code is generated automatically by pyR@TE, {} for the model : {}\n\"\"\"\n'.format(
                        date,model._Name)
        CodeStr += "from __future__ import division\nimport numpy as np\ndef beta_function_{}(t,y,Assumptions={{'two-loop': False,'diag': True}}):\n\t#The Assumption can be used to incorporate some switch in the calculation\n".format(model._Name)
        #We need an identity for each yukawa in principle because they could have different shape
        CodeStr += "\tkappa = 1./(16*np.pi**2)\n"
        #Simplify the keys of StrucYuk
        maps = {}
        for ikk,(key,val) in enumerate(StrucYuk.items()) :
                maps[key] = ''.join(reg.split('\\\(.*)',''.join(reg.split('{|}',key))))
                if max(val) != 1 :
                        CodeStr += "\tId{} = np.eye({})\n".format(ikk,max(val))
                else :
                        CodeStr += "\tId{} = 1\n".format(ikk)
        #We've got to get all the different RGEs and create variables
        NbRGEs = 0 
        #List of all the symbols
        ListSymbs = [el[1].g for el in model.GaugeGroups] + model.ListLbd + model.ListScM + model.ListTri
        for iel,el in enumerate(ListSymbs) :
                if len(reg.split('{(.*)}',str(el))) == 3 or len(reg.split('\\\(.*)',str(el))) == 3 :
                        ListSymbs[iel] = ''.join(reg.split('\\\(.*)',''.join(reg.split('{(.*)}',str(el)))))
                        maps[el] = ListSymbs[iel]       
                else :
                        maps[el] = el
        #First the scalar quatities
        NbRGEs += len(ListSymbs)
        #Now the non scalars quantities go through all of them and count the multiplicity due to the number of Generations
        if 'Yukawas' in model.Potential : 
                ToLoop = {}
                ToLoop.update(model.Potential['Yukawas'])
                if 'FermionMasses' in model.Potential : 
                        ToLoop.update(model.Potential['FermionMasses'])
                for key,val in ToLoop.items() : 
                        #Is it a squared matrix ? without the scalars
                        if len(StrucYuk[key]) == 3 :
                                if StrucYuk[key][-1] != 1 :
                                        exit("Index families on higgs are not taken into account")
                                else :
                                        StrucYuk[key] = StrucYuk[key][:-1]
                        #it can be scalar
                        if StrucYuk[key] == [1,1] :
                                NbYuk = 1
                                StrucYuk[key].append('scalar')
                        else :
                                NbYuk = StrucYuk[key][1]*StrucYuk[key][0]
                                StrucYuk[key].append('matrix')
                        NbRGEs += NbYuk
        #Count the rges that we have written 
        iwritten = 0
        #declare the matrices
        Mapping = {}
        ##########################################################################################
        #Aux function
        def declarematrices(llist,iwritten,Mapping,StrucYuk,CodeStr):
                for ii,val in llist:
                        if StrucYuk[ii][-1] == 'matrix' :
                                CodeStr += "\t{} = np.matrix(np.zeros(({},{})))\n".format(maps[ii],StrucYuk[ii][0],StrucYuk[ii][1])
                                Mapping[ii]=(iwritten,iwritten+StrucYuk[ii][0]*StrucYuk[ii][1]-1)
                                #fill up the newly constructed matrix
                                for j in range(StrucYuk[ii][0]):
                                        for k in range(StrucYuk[ii][1]):
                                                CodeStr += "\t{0}[{1},{2}] = y[{3}]\n".format(maps[ii],j,k,iwritten)
                                                iwritten += 1
                        elif StrucYuk[ii][-1] == 'scalar' : 
                                Mapping[ii]=(iwritten,iwritten)
                                #fill up the newly constructed matrix
                                CodeStr += "\t{} = y[{}]\n".format(maps[ii],iwritten)
                                iwritten += 1
                return Mapping,iwritten,CodeStr
        ##########################################################################################
        if 'Yukawas' in FinalExpr.keys() : 
                Mapping,iwritten,CodeStr = declarematrices(FinalExpr['Yukawas'].items(),iwritten,Mapping,StrucYuk,CodeStr)
        if 'FermionMass' in FinalExpr.keys() :
                Mapping,iwritten,CodeStr = declarematrices(FinalExpr['FermionMass'].items(),iwritten,Mapping,StrucYuk,CodeStr)
        #Do the mapping once and for all
        for el in ListSymbs:
                Mapping[el] = iwritten
                iwritten += 1
        CodeStr += "\t###############\n\t#beta functions\n\t###############\n"
        CodeStr += "\t#This is the mapping between rges and y[i]\n\t#{}\n".format(Mapping)
        temp1L=temp2L = ''
        for ii,val in FinalExpr.items() : 
                for label,expr in val.items() : 
                                if ii == 'Gauge-Couplings' :
                                        label = model.GetGroupFromName[label].g
                                #Get the 1loop contribution
                                if ii != 'Yukawas' and ii != 'FermionMass' : 
                                        temp1L += '\tb{} = ({})*kappa\n'.format(
                                                        maps[label],TranslateToNumerics((kappa*expr).expand().subs(1/pi**2,0),ListSymbs,maps[label],model,Mapping,isScalar=True))
                                        #Get the 2loop contribution
                                        temp2L += '\t\tb{0} = b{0} + ({1})*kappa**2\n'.format(
                                                        maps[label],TranslateToNumerics((kappa**2*expr).expand().subs(pi,0),ListSymbs,maps[label],model,Mapping,isScalar=True))
                                else :
                                        temp1L += '\tbeta{} = ({})*kappa\n'.format(
                                                        maps[label],TranslateToNumerics((kappa*expr).expand().subs(1/pi**2,0),ListSymbs,maps[label],model,Mapping,isScalar=False))
                                        #Get the 2loop contribution
                                        temp2L += '\t\tbeta{0} = beta{0} + ({1})*kappa**2\n'.format(
                                                        maps[label],TranslateToNumerics((kappa**2*expr).expand().subs(pi,0),ListSymbs,maps[label],model,Mapping,isScalar=False))
        #We regroup all the 1 loop and 2 loop contributions to avoid multiple if statements
        CodeStr += temp1L
        CodeStr += "\tif Assumptions['two-loop']:\n{}".format(temp2L)
        #Create the return and write the file
        def Sort(x) :
                if type(x[1]) == int : 
                        return x[1]
                else :
                        return x[1][0]
        SortMapping = sorted(Mapping.items(),key=Sort)
        Return,ReturnFull = "",""
        Tags,TagsFull = [],[]
        for el in SortMapping : 
                if type(el[1]) != int and StrucYuk[el[0]][-1] =='matrix': 
                        for ii in range(StrucYuk[el[0]][0]) : 
                                for jj in range(StrucYuk[el[0]][1]):
                                        if ii == jj :
                                                Return += 'beta{}[{},{}], '.format(maps[el[0]],ii,ii)
                                                Tags.append('\"beta{}{}\"'.format(maps[el[0]],ii,ii))
                                        else :
                                                #If it is a vector all component should be returned 
                                                if StrucYuk[el[0]][1] == 1  : #vector
                                                        Return += 'beta{}[{},{}], '.format(maps[el[0]],ii,jj)
                                                        Tags.append('\"beta{}[{},{}]\"'.format(maps[el[0]],ii,jj))
                                                else :
                                                        Return += '0,'
                                                        Tags.append('0')
                                        ReturnFull += 'beta{}[{},{}], '.format(maps[el[0]],ii,jj)
                                        TagsFull.append('\"beta{}[{},{}]\"'.format(maps[el[0]],ii,jj))
                elif type(el[1]) != int and StrucYuk[el[0]][-1] =='scalar':
                        Return += 'beta{}, '.format(maps[el[0]])
                        ReturnFull += 'beta{}, '.format(maps[el[0]])
                        TagsFull.append('\"beta{}\"'.format(maps[el[0]]))
                elif type(el[1]) == int : 
                        Return += 'b{}, '.format(el[0])
                        ReturnFull += 'b{}, '.format(el[0])
                        TagsFull.append('\"b{}\"'.format(el[0]))
                else :
                        exit("Error this is not implemented, contact the authors")
        CodeStr += "\tif Assumptions['diag'] :\n\t\treturn [{}]\n\telse :\n\t\treturn [{}]".format(Return[:-2],ReturnFull[:-2])#remove the last coma
        #transform Tag
        Tags= eval('[{}]'.format(','.join(Tags)))
        TagsFull= eval('[{}]'.format(','.join(TagsFull)))
        CodeStr +=      "\n\n\n\n#rge = RGE(beta_function_{},{},labels={}).\n".format(model._Name,NbRGEs,TagsFull)
        #Check that there are no points and that it ends with a .py
        if settings['ExportFile'] =='BetaFunction.py':
                settings['ExportFile'] = 'BetaFunction{}.py'.format(model._Name)
        spltName = settings['ExportFile'].split('.')
        if len(spltName) > 2 :
                Name = spltName[0]
                for ii in range(1,len(spltName)-1):
                        Name = Name + spltName[ii] 
                if spltName[-1] == 'py':
                        Name = Name+'.py'
                else :#we add the rest and add the dot py by force
                        Name = Name+spltName[-1] +'.py'
        elif len(spltName) == 1 :
                Name = spltName[0] + '.py'
        elif len(spltName) == 2 and spltName[-1] != 'py':
                Name = spltName[0] + spltName[1] +'.py'
        else :
                Name = spltName[0] + '.'+ spltName[1]
        f = open('{}/{}'.format(settings['Results'],Name),'w')
        f.write(CodeStr)
        f.close()
        #################################
        #End of the beta_function.py file
        #################################
        #Now let's create an example file in which we demonstrate how to use the above file in combination with the RGEClass.py
        #Header
        CodeStrI = '\"\"\"\nThis Code is generated automatically by pyR@TE, {} for the model : {}\nIt demonstrates how to use the beta_function file produced by pyR@TE using the -e option to solve the rges using the RGEsclass.py module\"\"\"\n'.format(
                        date,model._Name)
        CodeStrI += "#!/usr/bin/env python\nimport sys\n#This line should be modified if you move the result file\nsys.path.append('{}')\n#Same as above.The RGEsclass.py should be accessible\nsys.path.append('{}')\ntry :\n\tfrom RGEclass import *\nexcept ImportError :\n\texit('RGEclass.py not found.')\ntry :\n\timport scipy.interpolate as scp\nexcept ImportError :\n\texit('error while loading scipy')\nfrom {} import *\n".format(settings['Results'],'{}/Source/Output'.format(cwd),Name.split('.')[0]) 
        CodeStrI += "#In Order to use your result you have to create an instance of the RGE class following this\nrge = RGE(beta_function_{},{},labels={})\n#Note here that the labels are just used to label the results but should match the order in which the RGEs are solved.\n".format(model._Name,NbRGEs,TagsFull)
        CodeStrI += "\"\"\"\nOne can access the physical parameters as the top mass, the up quark yukawas, Z mass ...\nSince they are all attributes of the class object we have just created. The list is the following :\n\tPi = np.pi\n\talpha = 1./128.91\n\talphaS = 0.1184\n\tsinthetasq = 0.2316 #Sin squared\n\tMz = 91.1876 #in GeV\n\tMw = 80.399 # in GeV\n\tGf = 1.16637e-5 #in GeV ^-2\n\tCF = 4./3. #Casimir of the SU3 color\n\tNc = 3.\n\tQtop = 2./3.#top quark chargs\n\tvev = 246.22 # From (sqrt(2)Gf)^-1/2\n\tMt = 173.3 # the pole mass of the top\n\tU = [2.4e-3,1.27,165.3308]#mu,mc,mt (MSbar)mt=(Mt/(1+4alphas/3Pi))\n\tD = [4.75e-3,104e-3,4.19]#md,ms,mb\n\tL = [0.511e-3,105.66e-3,1.777]\n\tyU = [sqrt(2)/self.vev * el for el in self.U]\n\tyD = [sqrt(2)/self.vev * el for el in self.D]\n\tyL = [sqrt(2)/self.vev * el for el in self.L] \n\tg10 = sqrt(4*self.Pi * self.alpha)/sqrt(1-self.sinthetasq)#Initial value\n\tg20 = sqrt(4*self.Pi * self.alpha)/sqrt(self.sinthetasq)#Initial value\n\tg30 =  sqrt(4*self.Pi * self.alphaS)#Initial value\n\n\nTo modify one of them one would do :\n\"\"\"\nrge.Mt = 171.0\n"
        CodeStrI += "#Now one should specify the initial values for the {0} parameters of the rges\n########################\n#USER INPUT REQUIRED HERE\n########################\nY0 = [0.0]*{0}\n#assign this list to the initial values of the rges\nrge.Y0 = Y0".format(NbRGEs)
        CodeStrI += "\n\"\"\"the translation is the following : {}\"\"\"\n".format(Mapping)
        CodeStrI +="#Define the t values, using numpy: tmin,tmax,tstep.t=log_10(Scale in GeV)\ntmin,tmax,tstep = np.log(rge.Mz)/np.log(10),19,0.1\n"
        CodeStrI += "#Let's define the assumptions that can be passed through this dictionnary. If you want to add other than just two-loop you have to do the corresponding modification in {}.py\nAssumptions={{'two-loop':False,'diag': True}}\n".format(settings['ExportFile'])
        CodeStrI += "#Now we can finally solve the rges for the required values\nrge.solve_rges(tmin,tmax,tstep,Assumptions)\n"
        CodeStrI += "#Finally if one wants to access the result it can use the labels define above for an example see below. Note that the 't' entry is defined by default.\nfrom pylab import *\n#Some settings\nrc('axes',linewidth=2)\ngrid(True)\nfigprops = dict(figsize=(6., 5.))\nfigure(1,**figprops)\nadjustprops = dict(top=0.95, bottom=0.15, left=0.18, right=0.95)\nsubplots_adjust(**adjustprops)\nplot(rge.Sol['t'],rge.Sol['{3}'],label='{3} of {0} model')\nplot(rge.Sol['t'],rge.Sol['{1}'],label='{1} of {0} model')\nplot(rge.Sol['t'],rge.Sol['{2}'],label='{2} of {0} model')\n#Some more settings\nxlabel(r'$t$',fontsize =18)\nylabel('{1},{2} and {3} params of {0}')\nleg=gca().legend(ncol=3,loc=3,shadow=True,numpoints=1)\nshow()\n".format(model._Name,TagsFull[0],TagsFull[1],TagsFull[2])
        f = open('{}/{}'.format(settings['Results'],'SolveRGEs.py'),'w')
        f.write(CodeStrI)
        f.close()
        return '{}/{}'.format(settings['Results'],Name)



def ExportBetaToCpp(FileNumpyBeta,settings={}):
    """
    This is the export function of the beta function in Cpp function that can then be compiled into a c++ library that can then be Loaded with PyR@TE Solving wrapper.
    This is a script that translate the result of Numpy directly into Cpp compatible function using Template. The idea is that it should be possible to translate numpy functions into cpp without re-running the code (two-Loop G221 are taking for ever ~ 4 days 10 Gb of memory"""
    try :
        numpybeta = file(FileNumpyBeta,'r').readlines()
    except IOError : 
        print("Error could not read the BetaFunction.py file, skipping C++ export")
        pass
    names = [el for el in numpybeta if el.find('return') != -1 ][0]
    names = reg.split('[, \[\]]',names)
    names = [el for el in names if reg.search(el,"['0''1''2''3''4''5''6''7''8''9']") == None and el.find('return') ==-1 and el != '\n'] 
    #we need to keep the ordering while removing multiple occurences (caused a bug in the first version)
    duplicatenames = []
    for el in names : 
      if not(el in duplicatenames):
        duplicatenames.append(el)
    names = duplicatenames 
    #number of rges, whatch out each matrix is counted as one element
    nbrge = len(names) 
    MatrixNames = [el for el in numpybeta if el.find('np.matrix') != -1]
    Dim = [reg.search('.*np.zeros\(\((.),(.).*',el) for el in MatrixNames]
    MatrixNames = [reg.split('\t|=',el)[1] for el in MatrixNames]
    MatrixNames = [el.split(' ')[0] for el in MatrixNames]
    Init = [[reg.search(r'.*({})\[(.),(.)\][ ]*=[ ]*y\[(.{{1,2}})\].*'.format(el),line) for el in MatrixNames] for line in numpybeta]
    Init = [(el.group(1),el.group(2),el.group(3),el.group(4)) for el in sum(Init,[]) if el != None]
    Initindices = [el[-1] for el in Init]
    #At this point it might happen that scalar terms have been used directly to type in the equation like m_t in t' models
    #We have to look for them and declare them as doubles.
    HiddenScalars = [reg.search(r'([^=\t ]*)[ ]*=[ ]*y\[(.{1,2})\].*',line)for line in numpybeta]
    HiddenScalars = [(el.group(1),el.group(2)) for el in HiddenScalars if el != None]
    #keep only the ones that are not matrices 
    HiddenScalars = [el for el in HiddenScalars if not(el[1] in Initindices)]
    if len(Dim) != len(MatrixNames) : 
        print("error, while writing the interface to  c++, please report to the authors")
        pass
    else :
        Yuks = [(el,dim.group(1),dim.group(2)) for el,dim in zip(MatrixNames,Dim)]
        MatrixNBeta = ['beta{}'.format(mat) for mat in MatrixNames]
        #These ones have to be declared as double in C++
        for el in MatrixNBeta:
          names.remove(el)
        Eqs = []
        #Collect one and two loop beta functions for all the names
        for elem in names+MatrixNBeta : 
            Eqs.append([el for el in numpybeta if el.find(elem) != -1 ][:2])
        #Group together the one and two loop expressions
        OneL = '\n'.join([el[0] for el in Eqs])
        TwoL = '\n'.join([el[1] for el in Eqs])
        MatrixNBetaDeclare = ['double {} = 0.0;'.format(el) for el in names]
        #Get the info about all the matrices in the file
        #For each matrix declare it 
        Filling = []
        #MatrixDeclare.append("gsl_matrix * {} = gsl_matrix_alloc({},{});".format(el[0],el[1],el[2]) for el in Yuks]
        MatrixDeclare = ["mat {}({},{},fill::zeros);".format(el[0],el[1],el[2]) for el in Yuks]
        #Declare the hidden scalars as double
        MatrixNBetaDeclare += ["mat beta{}({},{},fill::zeros);".format(el[0],el[1],el[2]) for el in Yuks]
        #For each matrix create a Transpose a HC and a Conjugate copy
        MatrixDeclarehc =["mat {0}hc(conj({0}.t()));".format(el[0]) for el in Yuks]
        MatrixDeclareT =["mat {0}T({0}.st());".format(el[0]) for el in Yuks]
        MatrixDeclareC =["mat {0}C(conj({0}));".format(el[0]) for el in Yuks]
        #Do the transpose effectively
        for elem in Init : 
            #keep on filling it up We assume that the y[i] are real 
            Filling.append("{}({},{}) = y[{}];".format(elem[0],elem[1],elem[2],elem[3]))
        #Switch = "\tint two  0;\n\tif (TwoLoop){\n\t\tint two = 0;\n\t}\n\telse{\ntwo = 1;\n}\n"
        #Change the hidden scalars for their actual y[x] value
        for el in HiddenScalars :
            OneL = reg.sub(r'{}(?<!beta{})'.format(el[0],el[0]),'y[{}]'.format(el[1]),OneL)
            TwoL = reg.sub(r'{}(?<!beta{})'.format(el[0],el[0]),'y[{}]'.format(el[1]),TwoL)
        Tocpp = reg.sub(r'y\[(.{1,2})\]\*\*(.)','pow(y[\\1],\\2)',OneL),reg.sub(r'y\[(.{1,2})\]\*\*(.)','pow(y[\\1],\\2)',TwoL)
        #WARNING HERE ONLY REAL??? F. on November 20th I changed the def of hc and think it is tken into account --> check
        #For the fermion and scalar mass terms we get rid of the conjugation because it is incompatible with the double type of c
        ##############
        Tocpp = [reg.sub(r'np.transpose\(np.conjugate\((y\[.{1,2}\])\)\)','\\1',el) for el in Tocpp]
        Tocpp = [reg.sub(r'np.conjugate\(np.transpose\((y\[.{1,2}\])\)\)','\\1',el) for el in Tocpp]
        ################
        #"Hermitian"
        Tocpp = [reg.sub(r'np.transpose\(np.conjugate\(([^\)]*)\)\)','conj(\\1)',el) for el in Tocpp]
        Tocpp = [reg.sub(r'np.conjugate\(np.transpose\(([^\)]*)\)\)','\\1hc',el) for el in Tocpp]
        #Conjugate
        #################
        #For the fermion and scalar mass terms we get rid of the conjugation because it is incompatible with the double type of c
        Tocpp = [reg.sub(r'np.conjugate\((y\[.{1,2}\])\)','\\1',el)for el in Tocpp]
        #################
        Tocpp = [reg.sub(r'np.conjugate\(([^\)]*)\)','\\1C',el)for el in Tocpp]
        #Transpose
        Tocpp = [reg.sub(r'np.transpose\((y\[.{1,2}\])\)','\\1',el)for el in Tocpp]
        Tocpp = [reg.sub(r'np.transpose\(([^\)]*)\)','\\1T',el)for el in Tocpp]
        Tocpp = [reg.sub(r'np.trace','trace',el) for el in Tocpp]
        #Trace ^2 it does happen !
        Tocpp = [reg.sub(r'(.)\*\*\)(.*?)ecart',')\\1,)\\2ecart(wop',el[::-1])[::-1] for el in Tocpp]
        #kappa ^2
        Tocpp = [reg.sub(r'\*kappa\*\*(.)','*pow(kappa,\\1)',el) for el in Tocpp]
        Tocpp = [reg.sub(r'\n',';\n',el) for el in Tocpp]
        Tocpp = [reg.sub(r';\n;\n',';\n\n',el) for el in Tocpp]
        #Ensure the floats on fraction
        Tocpp = [reg.sub(r'([0-9]+)\*','\\1.0*',el) for el in Tocpp]
        Tocpp = [reg.sub(r'\/([0-9]+)','/\\1.0',el) for el in Tocpp]
        #Output of the beta function has to copied into f[]
        #Watch out in armadillo one must use regular brackets to access matrix elements !!
        TocppF=sum(sum([[['beta{}({},{})'.format(elem[0],i,j) for j in range(int(elem[1]))]for i in range(int(elem[2]))] for elem in Yuks],[]),[])
        TocppF += names
        LabelRGEs = ["{}".format(el) for el in TocppF]
        TocppF = ["f[{}] = {};".format(iel,el) for iel,el in enumerate(TocppF)]
        #We define a template funciton 
        template = """//This is automatically generated by pyR@TE. It is a c++ version of the calculated beta function that can be compiled into a .so Library that is then loaded via the BetaFunction class into python and solved using the regular solver of pyR@TE. See Makefile and BetaFunction class definition in RGEclass.py for more info.Using Armadillo library.\n\n\n#include<iostream>\n#define ARMA_DONT_USE_WRAPPER\n#include <armadillo>\nusing namespace std;\nusing namespace arma;\n\textern "C"{{\n\t\textern void beta_function( double t, double y[], double f[],bool TwoLoop);\n\t}}\n\n\nvoid  beta_function (double t, double y[], double f[], bool TwoLoop)\n{{\n\tdouble kappa = 1./(16*pow(3.14,2));\n\t//Matrix Declaration\n\t{}\n\t{}\n\t{}\n\t{}\n\t//beta Declaration\n\t{}\n\t{}\n\t//Beta functions one loop\n{}\n\n\t//Two loop beta functions\n\tif (TwoLoop){{\n\t{}\n}}\n\t//Copy the values in f for the return\n\t{}\n}}\nint main()\n{{\nreturn 0;\n}}""".format('\n\t'.join(MatrixDeclare),'\n\t'.join(Filling),'\n\t'.join(MatrixDeclarehc),'\n\t'.join(MatrixDeclareC),'\n\t'.join(MatrixDeclareT),'\n\t'.join(MatrixNBetaDeclare),Tocpp[0],Tocpp[1],'\n\t'.join(TocppF))
        OutputFileName = (FileNumpyBeta.split('.py')[0]+'.cpp').split('/')[-1]
        f = open('{}'.format(FileNumpyBeta.split('.py')[0]+'.cpp'),'w')
        f.write(template)
        #Now one has to write the Makefile for the corresponding files
        #Copied from nCTEQ to be tested 
        if platform.system() == 'Darwin':
            templateMakefile = """libbeta_{0}.dylib: {0}.o\n\tg++ -dynamiclib -o $@ $+ -llapack -lblas\n{0}.o: {0}.cpp\n\tg++ -fPIC -c -o $@ $<\nclean :\n\trm *.o\n\trm libbeta*.dylib\n\n#You might want to try with gcc instead.""".format(OutputFileName.split('.cpp')[0])
        else :
            templateMakefile = """libbeta_{0}.so: {0}.o\n\tg++ -shared -o $@ $+ -llapack -lblas\n{0}.o: {0}.cpp\n\tg++ -fPIC -c -o $@ $<\nclean :\n\trm *.o\n\trm libbeta*.so\n\n#You might want to try with gcc instead.""".format(OutputFileName.split('.cpp')[0])
        MakefileName = FileNumpyBeta.split('/')
        MakefileName[-1] = "Makefile"
        f = open('/'.join(MakefileName),'w')
        f.write(templateMakefile)
        f.close()
        #####
        #Write the example file
        #####
        template = """#!/usr/bin/env python\n\n#################\n#This file is automatically generated by PyR@TE and illustrate how to use the cpp beta function with the RGE class to solve the RGEs efficiently.\n#################\n\nimport sys\n\nsys.path.append('{0}/Source/Output/')\n\n#import the BetaFunction class\nfrom BetaFunctionCppclass import BetaFunctionCpp\n#import the RGEclass\nfrom RGEclass import *\n#Create an instance of this class\nbeta = BetaFunctionCpp('{1}',{2})\n#Create an instance of the RGE class\nrge = RGE(beta.betafunction,{2},labels={3})\n#set the Y0 vector\nY0 = [0.1]*{2}\n#set the scales t\nts = 0,40,40.0/100\n#initialize the Y0 attribute of the rge instance\nrge.Y0 = Y0\n#solve as in the python case\nrge.solve_rges(ts[0],ts[1],ts[2])\n#The results are stored in Sol dictionary\nprint(rge.Sol)""".format(cwd,'libbeta_{}.dylib'.format(OutputFileName.split('.cpp')[0]),len(TocppF),LabelRGEs)
        f = open('{}/SolveRGEsCpp.py'.format('/'.join(FileNumpyBeta.split('/')[:-1])),'w')
        f.write(template)
        f.close()



def ToNumpy(func,xpr,Mat=False):
        pos0 = xpr.find('{}('.format(func)) + len(func)
        lstring = list(xpr)
        # and now find the matching round bracket
        pos1 =  findclosingbracket(lstring,pos0+1)
        # change all the comas in between the two brackets
        while True : 
                tpfind = xpr[pos0:].find(',')
                pcoma = pos0 + tpfind if tpfind != -1 else len(lstring)-1
                if pcoma < pos1 : 
                        lstring[pcoma] = '*'
                        pos0 = pcoma+1 
                else :
                        break
        if not(Mat) :
                return ''.join(lstring).replace(func,'np.{}'.format(func[1:]),1)
        else :
                if func == 'SScalarProd' :
                        xpr = ''.join(lstring)
                        src = reg.search('SScalarProd\(([^\*]*)\*([^\)]*)\)',xpr)
                        if src != None :
                                xpr = xpr.replace('SScalarProd({}*{})'.format(src.group(1),src.group(2)),'np.trace(np.transpose({})*{})'.format(src.group(1),src.group(2)))
                        else :
                                exit("error while translating the to numpy notation")
                        return xpr
                else :
                        return ''.join(lstring).replace(func,'',1)


def TAdj(xpr) : 
        pos0 = xpr.find('Adj(')
        #replace only the first occurence
        xpr = xpr.replace('Adj','np.transpose(np.conjugate',1)
        lstring = list(xpr)
        pos1 = findclosingbracket(lstring,pos0+len('np.transpose(np.conjugate('))
        #add a bracket to close the additional one we introduced above
        lstring.insert(pos1+1,')')
        return ''.join(lstring)


def RemoveIndices(xpr,el,Yuk=True):
        #The flag is here to deal with the MatMul 
        if Yuk :
                #Removes the indices of a Yuk Symbol
                while el+'[' in xpr : 
                        pos0 = xpr.find(el+'[') + len(el)
                        lstring = list(xpr)
                        #remove the squar bracket
                        lstring[pos0] = ''
                        pos1 = xpr[pos0+1:].find(']')+pos0+1
                        #remove indices
                        lstring = lstring[:pos0] + lstring[pos1+1:]
                        xpr = ''.join(lstring)
        else :
                while 'MatMul(' in xpr : 
                        pos0 = xpr.find('MatMul(') + len('MatMul(')
                        lstring = list(xpr)
                        lstring[pos0] = ''
                        pos1 = findclosingbracket(lstring,pos0+1)
                        if lstring[pos1+1] == ',':
                                pos2 = xpr[pos1+1:].find(')') + pos1+1
                                lstring = lstring[:pos1+1]+lstring[pos2+1:]
                        elif lstring[pos1+1] == ')' and lstring[pos1] == ')':
                                lstring[pos1] = ''
                        else :
                                loggingDebug("Warning case not taken into account in RemoveIndices")
                                pass
                        xpr = ''.join(lstring)
                        xpr = xpr.replace('MatMul','matMul',1)
        return xpr


def TranslateToNumerics(expression,ListSymbs,label,model,Mapping,Matrices=True,isScalar=False) : 
        """Implements all the modification to the expression to be evaluated numerically and return a string"""
        #Use the ToMathematicaNotation to clean up the result a bit
        StrXpr = ToMathematicaNotation(expression,model,FlagSquare=False)#to avoid the '('->'['
        if not(isScalar) :
                StrXpr = StrXpr.replace('trace','ttrace').replace('conj','np.conjugate').replace('Tp','np.transpose').replace('ScalarProd','SScalarProd')
        else :
                StrXpr = StrXpr.replace('trace','ttrace').replace('conj','np.conjugate').replace('Tp','np.transpose').replace('ScalarProd','SScalarProd')
        while 'ttrace(' in StrXpr       :
                StrXpr = ToNumpy('ttrace',StrXpr)
        while 'Adj(' in StrXpr :
                StrXpr = TAdj(StrXpr)
        #Gauge-Couplings,Quartic,Trilinear,scalar masses
        for el in ListSymbs : 
                StrXpr = StrXpr.replace(str(el),'y[{}]'.format(Mapping[el]))
        #Now let's deal with the Yukawas and the indices 
        if Matrices :#result in matrix form
                #Remove all the indices
                for el in model.ListYukawa :
                        el =''.join(reg.split('\\\(.*)',''.join(reg.split('{(.*)}',el))))
                        StrXpr = RemoveIndices(StrXpr,el,Yuk=True)
                StrXpr = RemoveIndices(StrXpr,'MatMul',Yuk=False)
                while 'matMul(' in StrXpr :
                        StrXpr = ToNumpy('matMul',StrXpr,Mat=True)#The Mat is used to remove the tag of the function like ScalarProd -> ''
                while 'SScalarProd(' in StrXpr:
                        StrXpr = ToNumpy('SScalarProd',StrXpr,Mat=True)
        return StrXpr




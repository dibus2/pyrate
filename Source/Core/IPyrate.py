"""
F. February 17th 2015.
This is the first implementation of the interactive database querry for PyR@TE.
It allows one to check efficiently what is implemented what is not and how to specify which contraction it wants.

"""

try :
    import cmd
    import os
    import readline
    import rlcompleter
    from sympy import symbols,Symbol,Rational,sqrt,IndexedBase
    if 'libedit' in readline.__doc__:
        readline.parse_and_bind("bind ^I rl_complete")
    else :
        readline.parse_and_bind("tab: complete")
    import pickle
except : 
    exit("Error while loading one of the modules: `cmd, os, readline, rlcompleter, pickle`")
import pudb


class IdbquerryMissingArgument(Warning):
    """
    No arguments have been specified
    """
    def __init__(self,ttype=''):

        if ttype == '':
            print("Please, specify the required arguments!")
        elif ttype == 'gauge':
            print("Gauge group not implemented")
        elif ttype =='irrep':
            print("Irrep not implemented for the specified group")

class IdbquerryWrongFormat(Warning):
    """
    Wrong format of arguments
    """
    def __init__(self,istuple,onlygauge=False):
        if not(onlygauge): 
            if istuple :
                print("The second argument must be the dynkin label of an irrep represented by a list. E.g. [1,0]")
            else :
                print("The second argument must be an integer representing the dimension of the irrep")
        else :
            print("Only one argument required: the gauge group.")

class IdbquerryWrongFormatMultipleIrreps(Warning):
    """
    Wrong format of arguments
    """
    def __init__(self):
        print("The second argument must be the list of irreps to contract represented by their dynkin label E.g. [[1,0],[0,1]]")

class IdbquerryWrongFormatNbArg(Warning):
    """
    Wrong format of arguments
    """
    def __init__(self,nbarg):
        print("This function requires {} arguments".format(nbarg))


class IdbquerryUnkownContraction(Warning):
    """
    Wrong format of arguments
    """
    def __init__(self):
        print("This contraction is not in the databased.")



class Idbquerry(cmd.Cmd):

    prompt = 'Querry: '
    doc_header = 'doc_header'
    misc_header = 'misc_header'
    undoc_header = 'undoc_header'
    ruler = '-'
    last_output = ''
    a,b,c,d = map(IndexedBase,['a','b','c','d'])


    intro = "=========================================================================================================\n\t\t Interactive PyR@TE mode to querry the database.\n\t\t The database is built on Susyno v3. February 17th.\n\t\t IMPORTANT: All the results are therefore bounded by the implicit definitions of Susyno\n\t\t to which we refer the user for further information. \n\t\t Florian Lyonnet: flyonnet@smu.edu 2015\n=========================================================================================================\n"

    def cmdloop(self, intro=None):
        try :
            #welcome message
            print(self.intro)
            #load the database
            print("loading database of CGCs...")
            localdir = os.path.realpath(os.path.dirname(__file__))
            fdb = open(localdir+'/../GroupTheory/CGCs.pickle','r')
            self.db = pickle.load(fdb)
            fdb.close()
            self.extractinfo()
            print('db version: {}'.format(self.db['date']))
            print("done!")
            #Let's create the list for the autocompletion
        except :
            exit("Error while loading the database")
        return cmd.Cmd.cmdloop(self,'')

    def extractinfo(self):
        self.gaugegroupimplemented = self.db.keys()
        self.gaugegroupimplemented.remove('date')
        self.gaugegroupimplemented = [str(el) for el in self.gaugegroupimplemented]

    def do_prompt(self,line):
        "change the interactive prompt"
        self.prompt = line + ': '

    def do_shell(self,line):
        output = os.popen(line).read()
        print(output)
        self.last_output = output
    
    def do_Invariants(self,line):
        #returns the invariant for the given contraction and gauge group
        ls = line.split(' ')
        pudb.set_trace()
        try :
            if len(ls) == 1 and ls[0] == '':
                raise IdbquerryMissingArgument('gauge')
            elif len(ls) == 1 :
                if ls[1] in self.gaugegroupimplemented :
                    raise IdbquerryMissingArgument('irrep')
                else :
                    raise IdbquerryWrongFormatMultipleIrreps()
            elif len(ls)>3 :
                raise IdbquerryWrongFormatNbArg(2)
            else :
                #We have the proper number of arguments
                ls = (ls[0],tuple([tuple(el) for el in eval(ls[1])]))
                whichinvariant = len(ls[1])
                if type(ls[1][0]) != tuple :
                    raise IdbquerryWrongFormatMultipleIrreps
                if not(ls[0] in self.gaugegroupimplemented):
                    raise IdbquerryMissingArgument('gauge')
                if whichinvariant == 2 : 
                    if not(ls[1] in self.db[ls[0]]['Bilinear']) :
                        raise IdbquerryMissingArgument('irrep')
                    else :
                        print(self.construct_contraction(self.db[ls[0]]['Bilinear'][ls[1]],2))
                elif whichinvariant == 3:
                    if not(ls[1] in self.db[ls[0]]['Trilinear']):
                        raise IdbquerryMissingArgument('irrep')
                    else :
                        print(self.construct_contraction(self.db[ls[0]]['Trilinear'][ls[1]],3))
                elif whichinvariant == 4 :
                    if not(ls[1] in self.db[ls[0]]['Quartic']):
                        raise IdbquerryMissingArgument('irrep')
                    else :
                        print(self.construct_contraction(self.db[ls[0]]['Quartic'][ls[1]],4))
                else :
                    raise IdbquerryUnkownContraction()
        except (IdbquerryMissingArgument,IdbquerryWrongFormatMultipleIrreps,IdbquerryWrongFormatNbArg,IdbquerryUnkownContraction):
            pass
                                   
    def construct_contraction(self,contraction,whichinvariant):
        """
        contstruct a human readable output of the contraction a la susyno
        """
        if whichinvariant == 4 :
            contraction = sum([self.a[i]*self.b[j]*self.c[k]*self.d[l]*oo for (i,j,k,l,oo) in contraction])
        elif whichinvariant == 3 :
            contraction = sum([self.a[i]*self.b[j]*self.c[k]*oo for (i,j,k,oo) in contraction])
        elif whichinvariant == 2 :
            contraction = sum([self.a[i]*self.b[j]*oo for (i,j,oo) in contraction])
        else :
            raise IdbquerryUnkownContraction
        return contraction

    
    def do_Dynkin(self,line):
        return self.do_generic(line,'Dynkin')

    def complete_Dynkin(self,text,line,begidx,endidx):
        return self.complete_generic(text,line,begidx,endidx,'Dynkin')

    def do_Casimir(self,line):
        return self.do_generic(line,'Casimir')
    
    def complete_Casimir(self,text,line,begidx,endidx):
        return self.complete_generic(text,line,begidx,endidx,'Casimir')

    def do_DimR(self,line):
        return self.do_generic(line,'DynkinToDim')

    def complete_DimR(self,text,line,begidx,endidx):
        return self.complete_generic(text,line,begidx,endidx,'DynkinToDim')

    def do_GetDynk(self,line):
        return self.do_generic(line,'DimToDynkin',istuple=False)

    def complete_GetDynk(self,text,line,begidx,endidx):
        return self.complete_generic(text,line,begidx,endidx,'DimToDynkin',dim=True)
    
    def do_AdjR(self,line):
        return self.do_generic_onearg(line,'Adj')

    def complete_AdjR(self,text,line,begidx,endidx):
        return self.complete_onearg(text,line,begidx,endidx,'Adj')

    def do_FondR(self,line):
        return self.do_generic_onearg(line,'Fond')

    def complete_FondR(self,text,line,begidx,endidx):
        return self.complete_onearg(text,line,begidx,endidx,'Fond')
    

    def complete_generic(self,text,line,begidx,endidx,function,dim=False):
        #need to split it in two
        ls = line.split(' ')
        if len(ls) == 3 : 
            if not(ls[1] in self.db):
                completions = ls[1]+ls[2]
            else :
                if ls[2] == '':
                    if not(dim):
                        completions = [str(list(f)) for f in self.db[ls[1]][function].keys()]
                    else :
                        completions = [str(f) for f in self.db[ls[1]][function].keys()]
                else :
                    if not(dim):
                        completions = [str(list(f)) for f in self.db[ls[1]][function].keys() if str(list(f)).startswith(text)]
                    else :
                        completions = [str(f) for f in self.db[ls[1]][function].keys() if str(f).startswith(text)]
        elif len(ls) == 2 : 
            if ls[1] == '' :
                completions = self.gaugegroupimplemented
            else :
                completions = [f for f in self.gaugegroupimplemented if f.startswith(text)]
        else :
            completions = self.gaugegroupimplemented
        return completions


    def do_generic_onearg(self,line,function):
        args = line.split(' ')
        try :
            if len(args) == 1 and args[0] == '':
                raise IdbquerryMissingArgument('gauge')
            elif len(args) > 1 :
                raise IdbquerryWrongFormat(onlygauge=True)
            else :
                if args[0] in self.db:
                    print(self.db[args[0]][function])
                else :
                    raise IdbquerryMissingArgument('gauge')
        except (IdbquerryWrongFormat,IdbquerryMissingArgument):
            pass

    def complete_onearg(self,text,line,begidx,endidx,function):
        ls = line.split(' ')
        if len(ls) == 2 : 
            if ls[1] == '' :
                completions = self.gaugegroupimplemented
            else :
                completions = [f for f in self.gaugegroupimplemented if f.startswith(text)]
        else :
            completions = self.gaugegroupimplemented
        return completions

    def do_generic(self,line,function,istuple=True):
        args = line.split(' ')
        if istuple :
            tocheck = list
        else :
            tocheck = int
        try : 
            if not(type(eval(line.split(' ')[-1])) == tocheck) :
                raise IdbquerryWrongFormat(istuple)
        except IdbquerryWrongFormat:
            return
        try :
            if len(args) == 3 and args[0] =='SU2':
                args[1] = line[4:]
                args.pop(-1)
                #then most likely the user put a space in the True
            if len(args) == 2 :
                group,irrep = args
                if group and irrep :
                    if group in self.db:
                        if istuple:
                            if tuple(eval(irrep)) in self.db[group][function]:
                                print self.db[group][function][tuple(eval(irrep))]
                            else :
                                raise IdbquerryMissingArgument('irrep')
                        else :
                            if eval(irrep) in self.db[group][function]:
                                print list(self.db[group][function][eval(irrep)])
                            else :
                                raise IdbquerryMissingArgument('irrep')
                    else :
                        raise IdbquerryMissingArgument('gauge') 
            else : 
                raise IdbquerryMissingArgument()
        except IdbquerryMissingArgument:
            pass


    def do_EOF(self,line):
        return True

if __name__ == '__main__':
    
    Idbquerry().cmdloop()


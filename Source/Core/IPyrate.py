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
    def __init__(self,type=''):

        if type == '':
            print("Please, specify the required arguments!")
        elif type == 'gauge':
            print("Gauge group not implemented")
        elif type =='irrep':
            print("Irrep not implemented for the specified group")




class Idbquerry(cmd.Cmd):


    prompt = 'Querry: '
    doc_header = 'doc_header'
    misc_header = 'misc_header'
    undoc_header = 'undoc_header'
    ruler = '-'
    last_output = ''

    intro = "=========================================================================================================\n\t\t Interactive PyR@TE mode to querry the database.\n\t\t The database is built on Susyno v3. February 17th.\n\t\t Florian Lyonnet: flyonnet@smu.edu 2015\n=========================================================================================================\n"


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

    def do_Casimir(self,line):
        args = line.split(' ')
        try :
            if len(args) == 3 and args[0] =='SU2':
                args[1] = line[4:]
                args.pop(-1)
                #then most likely the user put a space in the True
            if len(args) == 2 :
                group,irrep = args
                if group and irrep :
                    if group in self.db:
                        if tuple(eval(irrep)) in self.db[group]['Casimir']:
                            print self.db[group]['Casimir'][tuple(eval(irrep))]
                        else :
                            raise IdbquerryMissingArgument('irrep')
                    else :
                        raise IdbquerryMissingArgument('gauge') 
            else : 
                raise IdbquerryMissingArgument()
        except IdbquerryMissingArgument:
            pass

    def complete_Casimir(self,text,line,begidx,endidx):
        #need to split it in two
        ls = line.split(' ')
        if len(ls) == 3 : 
            if not(ls[1] in self.db):
                completions = ls[1]+ls[2]
            else :
                if ls[2] == '':
                    completions = [str(list(f)) for f in self.db[ls[1]]['Casimir'].keys()]
                else :
                    completions = [str(list(f)) for f in self.db[ls[1]]['Casimir'].keys() if str(list(f)).startswith(text)]
        elif len(ls) == 2 : 
            if ls[1] == '' :
                completions = self.gaugegroupimplemented
            else :
                completions = [f for f in self.gaugegroupimplemented if f.startswith(text)]
        else :
            completions = self.gaugegroupimplemented
        return completions



        #if len(ls) >= 3 and ls[1] in self.gaugegroupimplemented:
        #    completions = [str(f) for f in self.db[ls[1]]['Casimir'].keys() if str(f).startswith(ls[2])]
        #elif len(ls) >=2 and ls[1] in self.gaugegroupimplemented:
        #    completions = [str(ls[1])+' '+ str(f) for f in self.db[ls[1]]['Casimir'].keys()]
        #else :
        #    completions = self.gaugegroupimplemented
        

#    def parseline(self,line):
#        parsed = line.split(' ')
#        return parsed[0],parsed[-1],line




    def do_EOF(self,line):
        return True

if __name__ == '__main__':
    Idbquerry().cmdloop()


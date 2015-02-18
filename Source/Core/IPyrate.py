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

    def do_prompt(self,line):
        "change the interactive prompt"
        self.prompt = line + ': '

    def do_shell(self,line):
        output = os.popen(line).read()
        print(output)
        self.last_output = output

    def do_Casimir(self,*args):
        try :
            if len(args[0]) == 2 :
                group,irrep = args[0]
                if group and irrep :
                    if group in self.db:
                        if eval(irrep) in self.db[group]['Casimir']:
                            print self.db[group]['Casimir'][eval(irrep)]
                        else :
                            raise IdbquerryMissingArgument('irrep')
                    else :
                        raise IdbquerryMissingArgument('gauge') 
            else : 
                raise IdbquerryMissingArgument()
        except IdbquerryMissingArgument:
            pass

#Does not work Fix it !

    def complete_Casimir(self,text,line,begidx,endidx):
        #need to split it in two
        pudb.set_trace()
        if len(line.split(' ')) >=2 and line.split(' ')[1] in self.gaugegroupimplemented:
            completions = [f for f in self.db[text]['DimToDynkin'].values() if f.startswith(text.split(' ')[1])]
        else :
            completions = self.gaugegroupimplemented
        return completions
#Always returns SU
#####


        

    def parseline(self,line):
        parsed = line.split(' ')
        return parsed[0],parsed[1:],line


    def do_EOF(self,line):
        return True

if __name__ == '__main__':
    Idbquerry().cmdloop()


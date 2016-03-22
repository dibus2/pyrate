"""
F. February 17th 2015.
This is the first implementation of the interactive database querry for PyR@TE.
It allows one to check efficiently what is implemented what is not and how to specify which contraction it wants.

"""

try:
    import cmd
    import os
    from sys import exit
    import sys
    sys.path.insert(0, '/Applications/HEPtools/sympy-0.7.6')
    import readline
    import rlcompleter
    from sympy import symbols, Symbol, Rational, sqrt, IndexedBase, Matrix
    import pickle
    import pudb
    import os
    sys.path.append('./../../../PyLie/git/')
    from PyLie import *
    import gzip
except:
    exit("Error while loading one of the modules: `cmd, os, readline, rlcompleter, pickle, gzip, PyLie`")


class IdbquerryMissingArgument(Warning):
    """
    No arguments have been specified
    """

    def __init__(self, ttype=''):

        if ttype == '':
            print("Please, specify the required arguments!")
        elif ttype == 'gauge':
            print("Gauge group not implemented")
        elif ttype == 'irrep':
            print("Irrep not implemented for the specified group")


class IdbquerryWrongFormat(Warning):
    """
    Wrong format of arguments
    """

    def __init__(self, istuple, onlygauge=False):
        if not (onlygauge):
            if istuple:
                print("The second argument must be the dynkin label of an irrep represented by a list. E.g. [1,0]")
            else:
                print("The second argument must be an integer representing the dimension of the irrep")
        else:
            print("Only one argument required: the gauge group.")


class IdbquerryWrongFormatMultipleIrreps(Warning):
    """
    Wrong format of arguments
    """

    def __init__(self):
        print(
            "The second argument must be the list of irreps to contract represented by their dynkin label E.g. [[1,0],[0,1]]")


class IdbquerryWrongFormatNbArg(Warning):
    """
    Wrong format of arguments
    """

    def __init__(self, nbarg):
        print("This function requires {} arguments".format(nbarg))


class IdbquerryUnkownContraction(Warning):
    """
    Wrong format of arguments
    """

    def __init__(self):
        print("This contraction is not in the databased.")


class IdbquerryInconsistentIrreps(Warning):
    """
    Wrong format of arguments
    """

    def __init__(self):
        print("The irrep should all be of the same gauge group i.e. have the same size.")


class Idbquerry(cmd.Cmd):
    prompt = 'Query: '
    doc_header = 'doc_header'
    misc_header = 'misc_header'
    undoc_header = 'undoc_header'
    ruler = '-'
    last_output = ''
    a, b, c, d = map(IndexedBase, ['a', 'b', 'c', 'd'])

    intro = "=========================================================================================================\n\t\t Interactive PyR@TE: mode to querry the database.\n\t\t The database is built on Susyno v3.4 November 6th.\n\t\t IMPORTANT: All the results are therefore bounded by the implicit definitions of Susyno\n\t\t to which we refer the user for further information. \n\t\t Florian Lyonnet: flyonnet@smu.edu 2015\n=========================================================================================================\n"

    def cmdloop(self, intro=None):
        try:
            # welcome message
            print(self.intro)
            # load the database
            print("loading database of CGCs...")
            self.localdir = os.path.realpath(os.path.dirname(__file__))
            with gzip.open(self.localdir + '/../GroupTheory/CGCs.pklz', 'rb') as fdb:
                self.db = pickle.load(fdb)
            fdb.close()
            self.extractinfo()
            print('db version: {}'.format(self.db['date']))
            print("done!")
            self.convert = {2: 'Bilinear', 3: 'Trilinear', 4: 'Quartic'}
            self._tosave = False
            # Let's create the list for the autocompletion
        except:
            exit("Error while loading the database")
        return cmd.Cmd.cmdloop(self, '')

    def do_exit(self, line):
        # store the db
        if self._tosave:
            self._update_db()
        exit("Exiting the interactive mode.")

    def _update_db(self):
        with gzip.open(self.localdir + '/../GroupTheory/CGCs.pklz', 'wb') as fdb:
            pickle.dump(self.db, fdb)
        fdb.close()

    def do_quit(self, line):
        if self._tosave:
            self._update_db()
        sys.exit("Exiting the interactive mode.")

    def extractinfo(self):
        self.gaugegroupimplemented = self.db.keys()
        if self.gaugegroupimplemented != []:
            self.gaugegroupimplemented.remove('date')
        self.gaugegroupimplemented = [str(el) for el in self.gaugegroupimplemented]

    def do_prompt(self, line):
        "change the interactive prompt"
        self.prompt = line + ': '

    def do_shell(self, line):
        output = os.popen(line).read()
        print(output)
        self.last_output = output

    def do_Invariants(self, line):
        # returns the invariant for the given contraction and gauge group
        ls = line.split(' ')
        if len(ls) > 2:
            ls = [ls[0], ''.join(ls[1:])]
        try:
            try:
                if len(ls) == 1 and ls[0] == '':
                    raise IdbquerryMissingArgument('gauge')
                elif len(ls) == 1:
                    if ls[1] in self.gaugegroupimplemented:
                        raise IdbquerryMissingArgument('irrep')
                    else:
                        raise IdbquerryWrongFormatMultipleIrreps()
                elif len(ls) > 3:
                    raise IdbquerryWrongFormatNbArg(2)
                else:
                    # We have the proper number of arguments
                    ls = (ls[0], tuple([tuple(el) for el in eval(ls[1])]))
                    whichinvariant = len(ls[1])
                    # Let's first the checks on the input and then validate wether it is in the db or not
                    if type(ls[1][0]) != tuple:
                        raise IdbquerryWrongFormatMultipleIrreps
                    if 'SU' in ls[0] and ls[0] != 'SU2':
                        if not (all([len(el) == int(ls[0][-1]) - 1 for el in ls[1]])):
                            raise IdbquerryInconsistentIrreps()
                    # Let's restrict it to SU(n) for now
                    if 'SU' not in ls[0]:
                        raise IdbquerryMissingArgument('gauge')
                    # Create the LieAlgebra object
                    try:
                        lie = LieAlgebra(CartanMatrix('SU', int(ls[0][-1])))
                    except ValueError:
                        exit("Error in creating the SU(n) gauge group. Possible groups are SU2,SU3,...")
                    if not self.is_in_db([ls[0], self.convert[whichinvariant], ls[1]]):
                        invs, conjs = self._split_inputs(ls[1])
                        res = lie.invariants(invs, conj=conjs, pyrate_normalization=True, returnTensor=True)
                        self._add_to_db([ls[0], self.convert[whichinvariant], ls[1]], res)
                        print(self.construct_contraction(res, whichinvariant))
                    else:
                        print(self.construct_contraction(self.db[ls[0]][self.convert[whichinvariant]][ls[1]],
                                                         whichinvariant))

            except (IdbquerryMissingArgument, IdbquerryWrongFormatMultipleIrreps, IdbquerryWrongFormatNbArg,
                    IdbquerryUnkownContraction, IdbquerryInconsistentIrreps):
                pass
        except SyntaxError:
            print("Could not interpret your querry, try again.")
            pass

    def construct_contraction(self, contraction, whichinvariant):
        """
        contstruct a human readable output of the contraction a la susyno
        """
        Contraction = ""
        if type(contraction[0]) != list:
            contraction = [contraction]
        for ii, ell in enumerate(contraction):
            if whichinvariant == 4:
                Contraction += "{} invariant:\n\t{}\n\n".format(ii + 1, sum(
                    [self.a[i] * self.b[j] * self.c[k] * self.d[l] * oo for (i, j, k, l, oo) in ell]))
            elif whichinvariant == 3:
                Contraction += "{} invariant:\n\t{}\n\n".format(ii + 1, sum(
                    [self.a[i] * self.b[j] * self.c[k] * oo for (i, j, k, oo) in ell]))
            elif whichinvariant == 2:
                Contraction += "{} invariant:\n\t{}\n\n".format(ii + 1,
                                                                sum([self.a[i] * self.b[j] * oo for (i, j, oo) in ell]))
            else:
                raise IdbquerryUnkownContraction
        return Contraction

    def do_Dynkin(self, line):
        return self.do_generic(line, 'dynkinIndex')

    def complete_Dynkin(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'Dynkin')

    def do_Casimir(self, line):
        return self.do_generic(line, 'casimir')

    def complete_Casimir(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'Casimir')

    def do_DimR(self, line):
        return self.do_generic(line, 'dimR')

    def complete_DimR(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'DynkinToDim')

    def do_GetDynk(self, line):
        return self.do_generic(line, 'DimToDynkin', istuple=False)

    def complete_GetDynk(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'DimToDynkin', dim=True)

    def do_AdjR(self, line):
        return self.do_generic_onearg(line, 'Adj')

    def complete_AdjR(self, text, line, begidx, endidx):
        return self.complete_onearg(text, line, begidx, endidx, 'Adj')

    def do_FondR(self, line):
        return self.do_generic_onearg(line, 'Fond')

    def complete_FondR(self, text, line, begidx, endidx):
        return self.complete_onearg(text, line, begidx, endidx, 'Fond')

    def complete_generic(self, text, line, begidx, endidx, function, dim=False):
        # need to split it in two
        ls = line.split(' ')
        if len(ls) == 3:
            if not (ls[1] in self.db):
                completions = ls[1] + ls[2]
            else:
                if ls[2] == '':
                    if not (dim):
                        completions = [str(list(f)) for f in self.db[ls[1]][function].keys()]
                    else:
                        completions = [str(f) for f in self.db[ls[1]][function].keys()]
                else:
                    if not (dim):
                        completions = [str(list(f)) for f in self.db[ls[1]][function].keys() if
                                       str(list(f)).startswith(text)]
                    else:
                        completions = [str(f) for f in self.db[ls[1]][function].keys() if str(f).startswith(text)]
        elif len(ls) == 2:
            if ls[1] == '':
                completions = self.gaugegroupimplemented
            else:
                completions = [f for f in self.gaugegroupimplemented if f.startswith(text)]
        else:
            completions = self.gaugegroupimplemented
        return completions

    def do_Matrices(self, line):
        mats = self.do_generic(line, 'Matrices', toreturn=True)
        print("Sparse matrices, showing non zero components:\n")
        print(mats)
        # print("{}".format("\n\n".join([str(Matrix(el)) for el in mats['mat']])))

    def complete_Matrices(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'Matrices')

    def do_generic_onearg(self, line, function):
        args = line.split(' ')
        try:
            if len(args) == 1 and args[0] == '':
                raise IdbquerryMissingArgument('gauge')
            elif len(args) > 1:
                raise IdbquerryWrongFormat(onlygauge=True)
            else:
                if args[0] in self.db:
                    print(list(self.db[args[0]][function]))
                else:
                    raise IdbquerryMissingArgument('gauge')
        except (IdbquerryWrongFormat, IdbquerryMissingArgument):
            pass

    def complete_onearg(self, text, line, begidx, endidx, function):
        ls = line.split(' ')
        if len(ls) == 2:
            if ls[1] == '':
                completions = self.gaugegroupimplemented
            else:
                completions = [f for f in self.gaugegroupimplemented if f.startswith(text)]
        else:
            completions = self.gaugegroupimplemented
        return completions

    def do_generic(self, line, function, istuple=True, toreturn=False):
        args = line.split(' ')
        if len(args) > 2:
            toremove = [(el + args[iel + 1], el, args[iel + 1], iel) for iel, el in enumerate(args) if el[-1] == ',']
            for ll in toremove:
                args.insert(ll[-1], ll[0])
                args.remove(ll[1])
                args.remove(ll[2])
            line = ' '.join(args)
        if istuple:
            tocheck = list
        else:
            tocheck = int
        try:
            if not (type(eval(line.split(' ')[-1])) == tocheck):
                raise IdbquerryWrongFormat(istuple)
        except IdbquerryWrongFormat:
            return
        try:
            if len(args) == 3 and args[0] == 'SU2':
                args[1] = line[4:]
                args.pop(-1)
                # then most likely the user put a space in the True
            if len(args) == 2:
                group, irrep = args
                irrep = eval(irrep)
                if istuple:
                    irrep = tuple(irrep)
                if group and irrep:
                    if not self.is_in_db([group, function, irrep]):
                        # use PyLie to calculate the corresponding result
                        assert 'SU' in group
                        if istuple:
                            if not len(irrep) == int(group[-1]) - 1:
                                raise IdbquerryInconsistentIrreps()
                            try:
                                lie = LieAlgebra(CartanMatrix("SU", int(group[-1])))
                            except ValueError:
                                exit("Error in creating the SU(n) gauge group. Possible groups are SU2,SU3,...")
                            # call the proper method from PyLie
                            res = eval('lie.{}({})'.format(function, irrep))
                            self._add_to_db([group, function, irrep], res)
                            if not toreturn:
                                print(res)
                            else:
                                return res
                        else:
                            # PyLie does not have a good implementation of the DimToDynkin therefore we will put it in the db directly
                            pudb.set_trace()
                            if function == 'DimToDynkin':
                                try:
                                    lie = LieAlgebra(CartanMatrix("SU", int(group[-1])))
                                except ValueError:
                                    exit("Error in creating the SU(n) gauge group. Possible groups are SU2,SU3,...")
                                # get all the irreps with dim <=
                                sign_irrep = sign(irrep)
                                res = [el for el in lie.repsUpToDimN(abs(irrep)) if lie.dimR(el) == abs(irrep)]
                                if res:
                                    if group != 'SU2':
                                        # TODO THERE IS A BUG IN THE SORTNIG OF REPSUPTODIMN FOR SU3 [1,1] e.g.
                                        if len(res) == 2:
                                            res = tuple(res[0]) if sign_irrep < 0 else tuple(res[1])
                                        else:
                                            res = tuple(res[0])  # For the adjoint for instance
                                    else:
                                        res = tuple(res[0]) if sign_irrep > 0 else (res[0][0], True)
                                else:
                                    res = []
                                self._add_to_db([group, function, irrep], res)
                                if not toreturn:
                                    print(res)
                                else:
                                    return res
                    else:
                        res = self.db[group][function][irrep]
                        if not toreturn:
                            print(res)
                        else:
                            return res
            else:
                raise IdbquerryMissingArgument()
        except IdbquerryMissingArgument:
            pass

    def do_EOF(self, line):
        return True

    # New for version 2 using directly PyLie
    def is_in_db(self, attributes):
        # Checks if the sequence of arguments are part of the db
        in_db = True
        for iel, el in enumerate(attributes):
            if iel == 0:
                if el not in self.db:
                    in_db = False
                    break
            elif iel == 1:
                if el not in self.db[attributes[0]]:
                    in_db = False
                    break
            elif iel == 2:
                if el not in self.db[attributes[0]][attributes[1]]:
                    in_db = False
                    break
            elif iel == 3:
                if el not in self.db[attributes[0]][attributes[1]][attributes[2]]:
                    in_db = False
                    break
        return in_db

    def _split_inputs(self, ls):
        invs, conjs = [], []
        for el in ls:
            check = [ell for ell in el if type(ell) == bool]
            if check != []:
                conjs.append(el[-1])
                invs.append(el[:-1])
            else:
                conjs.append(False)
                invs.append(el)
        return invs, conjs

    def _add_to_db(self, attributes, res):
        # add a given result to the database this way
        assert len(attributes) == 3
        self._tosave = True
        if attributes[0] in self.db:
            if attributes[1] in self.db[attributes[0]]:
                self.db[attributes[0]][attributes[1]][attributes[2]] = res
            else:
                self.db[attributes[0]][attributes[1]] = {attributes[2]: res}
        else:
            self.db[attributes[0]] = {attributes[1]: {attributes[2]: res}}






if __name__ == '__main__':
    Idbquerry().cmdloop()

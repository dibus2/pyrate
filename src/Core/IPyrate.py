"""
F. February 17th 2015.
This is the first implementation of the interactive database querry for PyR@TE.
It allows one to check efficiently what is implemented what is not and how to specify which contraction it wants.
F. April 2016: The interactive module has been merged with the main code of PyR@TE to access the group theory information
directly from PyLie.
"""

try:
    import cmd
    import os
    from sys import exit
    import sys
    import readline
    import rlcompleter
    import pickle
    import os
    from Logging import *
    import gzip
except:
    exit("Error while loading one of the modules: `cmd, os, readline, rlcompleter, pickle, gzip, PyLie`")
    sys.path.append(os.path.relpath(os.path.dirname(__file__)+'/../pylie'))
    from PyLie import *
try:
    Version = False
    from sympy import symbols, Symbol, Rational, sqrt, IndexedBase, Matrix, \
        Wild, Symbol, Function, symbols, pi, Rational, zeros, I, sqrt, eye, MatrixSymbol, \
        KroneckerDelta, flatten, pprint, IndexedBase, Idx, Integer, Add, Mul, Indexed, Sum, conjugate, adjoint, \
        __version__, Mod

    from sympy import MutableMatrix as tMM
    from sympy.physics.secondquant import evaluate_deltas
    #######
    # Shadowing of the piecewise mechanism of sympy that makes the code way too slow
    import sympy.concrete.expr_with_limits
    import sympy.concrete.delta
    from sympy.parsing.sympy_parser import parse_expr

    def mypiecewise_fold(exp):
        return exp
    sympy.concrete.expr_with_limits.piecewise_fold = mypiecewise_fold
    sympy.concrete.delta.piecewise_fold = mypiecewise_fold
    #######
    if __version__ != '1.0':
        Version = True
        raise ImportError
except ImportError:
    if Version:
        loggingCritical(
            "\tsympy version incompatible : {}, please get 1.0.0".format(__version__),
            verbose=True)
        exit()
    else:
        loggingCritical("\tError while loading sympy. Check the manual for required modules.", verbose=True)
        exit()

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

    intro =  """
                                              ___      _    _
                                             | _ \_  _| |  (_)___
                                             |  _/ || | |__| / -_)
                                             |_|  \_, |____|_\___|
                                                  |__/
        \n\n\n"""
    intro += "\t==============================================================================================\n\t\tPyLie: Python and Lie Algebra. PyLie has been written specifically for PyR@TE.\n\t\tIt is a python re-write of the relevant Mathematica routines of Susyno v3.4.\n\t\tFlorian Lyonnet: flyonnet@smu.edu 2015\n\t==============================================================================================\n"

    def __init__(self, noprint=False):
        cmd.Cmd.__init__(self)
        try:
            self.noprint = noprint
            # welcome message
            if not noprint:
                print(self.intro)
                # load the database
                print("loading database of CGCs...")
            self.localdir = os.path.realpath(os.path.dirname(__file__))
            with gzip.open(self.localdir + '/../GroupTheory/CGCs.pklz', 'rb') as fdb:
                self.db = pickle.load(fdb)
            fdb.close()
            self.extractinfo()
            if not noprint:
                print('db version: {}'.format(self.db['date']))
                print("done!")
            self.convert = {2: 'Bilinear', 3: 'Trilinear', 4: 'Quartic'}
            self._tosave = False
            # Let's create the list for the autocompletion
        except:
            exit("Error while loading the database")

    def cmdloop(self, intro=None):
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

    def do_Invariants(self, line, tensor=False):
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
                    # Let's restrict it to SU(n) and SO(N)
                    if ls[0][:2] not in ['SU', 'SO']:
                        raise IdbquerryMissingArgument('gauge')
                    # Create the LieAlgebra object
                    lie = self._declare_algebra(ls[0])
                    if not self.is_in_db([ls[0], self.convert[whichinvariant], ls[1]]):
                        invs, conjs = self._split_inputs(ls[1])
                        res = lie.invariants(invs, conj=conjs, pyrate_normalization=True, returnTensor=True)
                        self._add_to_db([ls[0], self.convert[whichinvariant], ls[1]], res)
                        if not self.noprint:
                            print(self.construct_contraction(res, whichinvariant))
                        else:
                            if not tensor:
                                return self.construct_contraction(res, whichinvariant)
                            else:
                                return res
                    else:
                        if not self.noprint:
                            print(self.construct_contraction(self.db[ls[0]][self.convert[whichinvariant]][ls[1]],
                                                         whichinvariant))
                        else:
                            if not tensor:
                                return self.construct_contraction(self.db[ls[0]][self.convert[whichinvariant]][ls[1]],
                                                         whichinvariant)
                            else:
                                return self.db[ls[0]][self.convert[whichinvariant]][ls[1]]

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
        if contraction == []:
            return "[]"
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
        return self.do_generic(self._filter_su2_irrep(line), 'dynkinIndex')

    def complete_Dynkin(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'Dynkin')

    def do_Casimir(self, line):
        # I need to protect these against calls with conjugated SU2 irrep e.g. [2,True]
        return self.do_generic(self._filter_su2_irrep(line), 'casimir')

    def complete_Casimir(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'Casimir')

    def do_DimR(self, line):
        return self.do_generic(self._filter_su2_irrep(line), 'dimR')

    def complete_DimR(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'DynkinToDim')

    def do_GetDynk(self, line):
        return self.do_generic(line, 'DimToDynkin', istuple=False)

    def complete_GetDynk(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'DimToDynkin', dim=True)

    def do_AdjR(self, line):
        return self.do_generic_onearg(line, 'adjoint')

    def complete_AdjR(self, text, line, begidx, endidx):
        return self.complete_onearg(text, line, begidx, endidx, 'Adj')

    def do_FondR(self, line):
        return self.do_generic_onearg(line, 'fond')

    def do_DimAdj(self, line):
        return self.do_generic_onearg(line, 'dimAdj')

    def do_ReduceProduct(self, line):
        self.noprint = True
        res = self.do_generic(line, 'reduceRepProduct')
        if res:
            res = [(el, self.do_generic(line.split(' ')[0]+' '+str(el[0]), 'dimR')) for el in res]
            res = sorted(res, key=lambda x: x[-1])
            print(res)
        self.noprint = False

    def complete_ReduceProduct(self, text, line, begidx, endidx):
        return self.complete_generic(text, line, begidx, endidx, 'reduceRepProduct')

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

    def do_Matrices(self, line, HB=False):
        if 'SU2' in line:
            line = self._filter_su2_irrep(line, domatrices=True)
        if not self.noprint:
            print("Sparse matrices, showing non zero components:\n")
        if type(line) == tuple and len(line) == 2:  # allows using conjugated irreps with SU2
            return self.do_generic(line[0], 'repMatrices', conjugate_matrix=line[-1], HB=HB)
        else:
            return self.do_generic(line, 'repMatrices', HB=HB)

    def do_Struc(self, line):
        return self.do_generic_onearg(line, '_get_struct()')

    def complete_Struc(self, text, line, begidx, endidx):
        return self.complete_onearg(text, line, begidx, endidx, '_get_struct()')

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
                if not self.is_in_db([args[0], function]):
                    lie = self._declare_algebra(args[0])
                    if function == '_get_struct()' or function == 'dimAdj':
                        res = eval('lie.{}'.format(function))
                    else:
                        res = eval('tuple(lie.{}.tolist()[0])'.format(function))
                    self._add_to_db([args[0], function], res)
                    if not self.noprint:
                        print(res)
                    else:
                        return res
                else:
                    if function in ['_get_struct()', 'dimAdj']:
                        if not self.noprint:
                            print(self.db[args[0]][function])
                        else:
                            return self.db[args[0]][function]
                    else:
                        if not self.noprint:
                            print(list(self.db[args[0]][function]))
                        else:
                            return list(self.db[args[0]][function])
        except (IdbquerryWrongFormat, IdbquerryMissingArgument):
            pass

    def do_PermutationSymmetryOfInvs(self, line):
        # Returns the permutation symmetry of the invariants under the form: (group1,group2), ((irrepSn1),mul), ((irrepSn2),mul))
        return self.do_generic(line, 'permutationSymmetryOfInvariants')

    def do_SnIrrepDim(self, line):
        args = line.split('SnDimIrrep')
        try:
            if len(args) != 1:
                raise IdbquerryMissingArgument()
            else:
                args = ''.join(args).replace(' ','')
                sn = Sn()
                try:
                    res = sn.snIrrepDim(eval(args))
                    print(res)
                except SyntaxError:
                    raise IdbquerryWrongFormat(args)
        except (IdbquerryWrongFormat, IdbquerryMissingArgument):
            return

    def do_RepsUpTo(self, line):
        args = line.split(' ')
        if len(args) != 2:
            raise IdbquerryMissingArgument()
        group, dim = args
        lie = self._declare_algebra(group)
        # call the proper method from PyLie
        res = lie.repsUpToDimN(int(dim))
        res = [(lie.dimR(el), el) for el in res]
        print(res)


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

    def do_generic(self, line, function, istuple=True, conjugate_matrix=None, HB=False):
        args = line.split(' ')
        if len(args) > 2 and function in ['dimR', 'casimir', 'dynkinIndex']:
            args = [args[0]] + [''.join(args[1:])]
        if len(args) > 2:
            toremove = [(el + args[iel + 1], el, args[iel + 1], iel) for iel, el in enumerate(args) if el[-1] == ',']
            for ll in toremove:
                args.insert(ll[-1], ll[0])
                args.remove(ll[1])
                args.remove(ll[2])
            line = ' '.join(args)
        elif len(args) == 1:
            try:
                raise IdbquerryMissingArgument()
            except IdbquerryMissingArgument:
                return None
        if istuple:
            tocheck = list
        else:
            tocheck = int
        try:
            try:
                if not (type(eval(args[-1])) == tocheck):
                    raise IdbquerryWrongFormat(istuple)
            except SyntaxError:
                raise IdbquerryWrongFormat(args[-1])
        except IdbquerryWrongFormat:
            return None
        try:
            if len(args) == 3 and args[0] == 'SU2':
                args[1] = line[4:]
                args.pop(-1)
                # then most likely the user put a space in the True
            if len(args) == 2:
                group, irrep = args
                irrep = eval(irrep)
                if istuple:
                    if type(irrep[0]) == list:
                        irrep = tuple([tuple(elem) for elem in irrep])
                    else:
                        irrep = tuple(irrep)
                if group and irrep:
                    if conjugate_matrix is not None:
                        irrep_tag = (irrep[0], conjugate_matrix)
                    else:
                        irrep_tag = irrep
                    if not self.is_in_db([group, function, irrep_tag]):
                        # use PyLie to calculate the corresponding result
                        if istuple:
                            lie = self._declare_algebra(group, irrep=irrep)
                            # call the proper method from PyLie
                            res = eval('lie.{}({})'.format(function, list(irrep)))
                            if function == 'repMatrices':
                                if conjugate_matrix is not None:
                                    if conjugate_matrix:
                                        res = [-el.conjugate() for el in res]
                                    irrep = (irrep[0], conjugate_matrix)
                                res = self._build_sparse_matrices(res)
                            self._add_to_db([group, function, irrep], res)
                            if not self.noprint:
                                if not function == 'repMatrices':
                                    print(res)
                                else:
                                    if not HB:
                                        print(res[0])
                                    else:
                                        print(res[1])
                            else:
                                if not function == 'repMatrices':
                                    return res
                                else:
                                    if not HB:
                                        return res[0]
                                    else:
                                        return res[1]
                        else:
                            # PyLie does not have a good implementation of the DimToDynkin therefore we will put it in the db directly
                            if function == 'DimToDynkin':
                                lie = self._declare_algebra(group)
                                # get all the irreps with dim <=
                                sign_irrep = sign(irrep)
                                res = sorted([el for el in lie.repsUpToDimN(abs(irrep)) if lie.dimR(el) == abs(irrep)])
                                if res:
                                    if group != 'SU2':
                                        if len(res) == 2:
                                            res = tuple(res[0]) if sign_irrep < 0 else tuple(res[1])
                                        elif len(res) == 1:
                                            res = tuple(res[0])
                                        else:
                                            pass  # multiple answers
                                    else:
                                        res = tuple(res[0]) if sign_irrep > 0 else (res[0][0], True)
                                else:
                                    res = []
                                self._add_to_db([group, function, irrep], res)
                                if not self.noprint:
                                    print(res)
                                else:
                                    return res
                    else:
                        if function == 'repMatrices':
                            if len(irrep) == 1:
                                irrep = (irrep[0], conjugate_matrix)
                            if not HB:
                                res = self.db[group][function][irrep]
                            else:
                                res = self.db[group]['HBmat'][irrep]
                        else:
                            res = self.db[group][function][irrep]
                        if not self.noprint:
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

    def _declare_algebra(self, group, irrep=[]):
        """
        creates the algebra depending on the info
        :param infos:
        :return:
        """
        if 'SU' in group:
            if irrep:
                if type(irrep[0]) == tuple or type(irrep[0]) == list:
                    for ell in irrep:
                        if not len(ell) == int(group[2:]) -1:
                            raise IdbquerryInconsistentIrreps()
                elif not len(irrep) == int(group[2:]) - 1:
                    raise IdbquerryInconsistentIrreps()
            try:
                lie = LieAlgebra(CartanMatrix("SU", int(group[2:])))
            except ValueError:
                    exit("Error in creating the SU(n) gauge group. Possible groups are SU2,SU3,...")
        elif 'SO' in group:
            try:
                lie = LieAlgebra(CartanMatrix("SO", int(group[2:])))
            except ValueError:
                exit("Error in creating the SO(n) gauge group.")
        else:
            exit("Exit group not implemented yet: `{}`".format(group[:2]))
        return lie

    def _add_to_db(self, attributes, res):
        # add a given result to the database this way
        self._tosave = True
        if len(attributes) == 3:
            if attributes[0] in self.db:
                if attributes[1] in self.db[attributes[0]]:
                    if attributes[1] == 'repMatrices':
                        self.db[attributes[0]][attributes[1]][attributes[2]] = res[0]
                        self.db[attributes[0]]['HBmat'][attributes[2]] = res[1]
                    else:
                        self.db[attributes[0]][attributes[1]][attributes[2]] = res
                else:
                    if attributes[1] == 'repMatrices':
                        self.db[attributes[0]][attributes[1]] = {attributes[2]: res[0]}
                        self.db[attributes[0]]['HBmat'] = {attributes[2]: res[1]}
                    else:
                        self.db[attributes[0]][attributes[1]] = {attributes[2]: res}
            else:
                if attributes[1] == 'repMatrices':
                    self.db[attributes[0]] = {attributes[1]: {attributes[2]: res[0]}}
                    self.db[attributes[0]]['HBmat'] = {attributes[2]: res[1]}
                else:
                    self.db[attributes[0]] = {attributes[1]: {attributes[2]: res}}
        elif len(attributes) == 2:
            if not attributes[0] in self.db:
                self.db[attributes[0]] = {attributes[1]: res}
            else:
                self.db[attributes[0]][attributes[1]] = res
        else:
            exit("Error in updating the db")

    def _build_sparse_matrices(self, res):
        sparse = []
        tp = {}
        shp = res[0].shape
        tp[(0, 0)] = {'mat': [], 'shape': shp}
        tp[(0, 1)] = {'mat': [], 'shape': shp}
        tp[(1, 0)] = {'mat': [], 'shape': shp}
        tp[(1, 1)] = {'mat': [], 'shape': shp}
        for mat in res:
            # Construct the matrix
            tp[(0, 0)]['mat'].append(Rational(1, 2) * (Matrix(mat) - Matrix(mat).conjugate()))
            tp[(0, 1)]['mat'].append(Rational(1, 2) * I * (Matrix(mat) + Matrix(mat).conjugate()))
            tp[(1, 0)]['mat'].append(-Rational(1, 2) * I * (Matrix(mat) + Matrix(mat).conjugate()))
            tp[(1, 1)]['mat'].append(Rational(1, 2) * (Matrix(mat) - Matrix(mat).conjugate()))

            for akey, aval in tp.items():
                tp[akey]['mat'][-1] = tuple(aval['mat'][-1])
        # reorganize the matrices we want in [(0,0)] all the matrices
        tp[(0, 0)]['mat'], tp[(0, 1)]['mat'], tp[(1, 0)]['mat'], tp[(1, 1)]['mat'] = tuple(tp[(0, 0)]['mat']), tuple(
            tp[(0, 1)]['mat']), tuple(tp[(1, 0)]['mat']), tuple(tp[(1, 1)]['mat'])
        # Now let's transform them to sparse matrices
        for struc, vstruc in tp.items():
            AllMat = []
            for matrix in vstruc['mat']:
                tempMat = []
                tpmatrix = np.array(matrix).reshape((vstruc['shape'][0], vstruc['shape'][1]))
                for line in range(vstruc['shape'][0]):
                    for col in range(vstruc['shape'][1]):
                        if tpmatrix[line][col] != 0:
                            tempMat.append((line, col, tpmatrix[line][col]))
                AllMat.append(tuple(tempMat))
            tp[struc] = cp.deepcopy(tuple(AllMat))

        for matrix in res:
            tempMat = []
            for line in range(matrix.shape[0]):
                for col in range(matrix.shape[1]):
                    if matrix[line, col] != 0:
                        tempMat.append((line, col, matrix[line, col]))
            sparse.append(tuple(tempMat))
        res = cp.deepcopy(tuple(sparse)), tp
        return res

    def toline(self, attributes):
        temp = [str(el) for el in attributes]
        if len(attributes) >=2:
            return ' '.join([temp[0],''.join(temp[1:])]).replace('(', '[').replace(')', ']').replace(', ',',')
        else:
            return ' '.join(temp).replace('(', '[').replace(')', ']')

    def _filter_su2_irrep(self, line, domatrices=False):
        if 'SU2' in line:
            ls = line.split(' ')
            if 'True' in ls[-1]:
                ls[-1] = ls[-1].replace('True', '')
                line = ' '.join(ls)
                if domatrices:
                    return line, True
                else:
                    return line
            if 'False' in ls[-1]:
                ls[-1] = ls[-1].replace('False', '')
                line = ' '.join(ls)
                if domatrices:
                    return line, False
                else:
                    return line
            else:
                if domatrices:
                    return line, False
                else:
                    return line
        else:
            return line


if __name__ == '__main__':
    Idbquerry().cmdloop()

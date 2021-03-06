#!/usr/bin/env python
from Contraction import *

""" In this file we define all tdhe symbols and conS*nstants of the groups that are used i.e. SU(N) U(1)"""
m, n, o, p, q, r, s, t, u, v = map(Wild, ['m', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v'])


class SUn(object):
    """Dynamical SU(N) class"""

    def __init__(self, N, name, idb):
        self.idb = idb  # This is the interactive db access object
        self.name = name
        self._absname = 'SU{}'.format(N)
        self.Algebra = 'A'
        self.N = N
        # Repr of the algebra or the fundamental
        self.Fond = self.idb.do_FondR(self.idb.toline([self._absname]))
        self.Adj = self.idb.do_AdjR(self.idb.toline([self._absname]))
        # dim of the adjoint
        self.dimAdj = self.idb.do_DimAdj(self._absname)
        self.g = symbols('g_{}'.format(self.name))
        self.singlet = '1'
        self.Dynksinglet = tuple((self.N - 1) * [0])
        self.fabc = self.idb.do_Struc(self.idb.toline([self._absname]))
        # transformation
        # Convert the  matrices into tuples that can be passed to the symbolic functions
        self.fabc = (tuple(flatten(self.fabc)), self.fabc.shape)
        self.U = False

    def dimR(self, irrep):
        return self.idb.do_DimR(self.idb.toline([self._absname, irrep]))

    def __repr__(self):
        """Change the representation of SU(N)"""
        return "Instance of SU({}) group".format(self.N)

    def MatRepr(self, ferm1, ferm2, HB=False):
        """Get the representation of ferm1,ferm2 """
        # The assertion here is that ferm1 and ferm2 are two different components of the same particle
        if not (HB):
            assert ferm1 == ferm2
            mat = self.idb.do_Matrices(
                self.idb._filter_su2_irrep(
                    self.idb.toline([self._absname, ferm1]), domatrices=True
                ), HB=False
            )

        else:
            assert ferm1[1] == ferm2[1]  # i.e. same representation
            # read the sequence i.e. PiPi PiSig SigPi
            mat = self.idb.do_Matrices(
                self.idb._filter_su2_irrep(
                    self.idb.toline([self._absname, ferm1[1]]), domatrices=True
                ), HB=True
            )[(ferm1[0], ferm2[0])]
        return mat

    def C2(self, irrep):
        """
        Casimir for the different irepp of SU(N)
        G stands for the adjoint given by the structure constants.
        """
        if irrep == self.Dynksinglet:
            return 0
        elif irrep != 'G':
            return self.idb.do_Casimir(self.idb.toline([self._absname, irrep]))
        elif irrep == 'G':
            return self.idb.do_Casimir(self.idb.toline([self._absname, self.Adj]))
        else:
            loggingCritical("Error: no such irrep {}".format(irrep), verbose=True)

    def S2(self, irrep):
        """Dynkin index for the SU(N) group. For the adjoint we have S2(G)=C2(G) see Langacker p110 for instance"""
        # translate the Dynkin notation
        if irrep == self.Dynksinglet:
            return 0
        elif irrep != 'G':
            return self.idb.do_Dynkin(self.idb.toline([self._absname, irrep]))
        elif irrep == 'G':
            return self.idb.do_Dynkin(self.idb.toline([self._absname, self.Adj]))
        else:
            loggingCritical("Error: no such irrep {}".format(irrep), verbose=True)


# New as of May 2016 Implementation of SO(2n) groups


class SOn(SUn):
    def __init__(self, N, name, idb):
        self.idb = idb
        self.name = name
        # The type of Lie algebra needs to be determined
        self.N = N
        if self.N % 2 == 0:
            self.Algebra = 'D'  # SO(2n)
        else:
            self.Algebra = 'B'  # SO(2n+1)
        self._absname = 'SO{}'.format(N)
        self.n = self.N / 2  # small n
        # dim of the adjoint
        self.Dynksinglet = tuple((self.n) * [0])
        # Repr of the algebra or the fundamental
        self.Fond = self.idb.do_FondR(self.idb.toline([self._absname]))
        self.Adj = self.idb.do_AdjR(self.idb.toline([self._absname]))
        self.dimAdj = self.idb.do_DimAdj(self._absname)
        self.fabc = self.idb.do_Struc(self.idb.toline([self._absname]))
        self.g = symbols('g_{}'.format(self.name))
        self.singlet = '1'
        self.U = False

    def __repr__(self):
        """Change the representation of SO(N)"""
        return "Instance of SO({}) group".format(self.N)


class U1(object):
    """
    Defines the U1 group class.
    The Casimir operator of the U1 Group is given by the hypercharge squared over 4. Y**2/4
    """

    def __init__(self, name):
        self.dimAdj = 1
        self.name = name
        if name != "U1":
            self.g = symbols('g_{}'.format(self.name))
        else:
            self.g = Symbol('g1')
        self._absname = 'U1'
        self.Algebra = ''
        self.singlet = '0'
        self.Dynksinglet = 0
        self.U = True
        self.fabcd = 0  # no structure constants

    def __repr__(self):
        """Change the representation of U1 instances."""
        return "Instance of the U1 group."""

    def Structure(self, a, b, c, d):
        return "WARNING U1 factors are abelian i.e. no structure constants"

    def MatRepr(self, f1, f2, HB=False):
        """Representation matrices of U1 group. This is for fermions. For Scalars see HB forms."""
        if not (HB):
            f1Ind, f2Ind = f1[1:], f2[1:]
            Indices = [[f1Ind[i], f2Ind[i]] for i in range(len(f1Ind))]
            return f1[0] * functools.reduce(operator.mul, [KroneckerDelta(i, j) for (i, j) in Indices], 1)
        else:
            assert f1[1] == f2[1]  # hypercharge
            Sequence = (f1[0], f2[0])
            f1Ind, f2Ind = f1[2:], f2[2:]
            # Sanity check
            if not (len(f1Ind) == len(f2Ind)):
                loggingCritical("Error, inconsistent value detected, contact the author")
                exit()
            # Indices = [[f1Ind[i],f2Ind[i]] for i in range(len(f1Ind))]
            Indices = zip(f1Ind, f2Ind)
            if Sequence == (0, 0):
                return Rational(1, 2) * (f1[1] - f1[1].conjugate()) * functools.reduce(operator.mul,
                                                                                       [KroneckerDelta(i, j) for (i, j)
                                                                                        in Indices], 1)
            elif Sequence == (0, 1):
                return Rational(1, 2) * I * (f1[1] + f1[1].conjugate()) * functools.reduce(operator.mul,
                                                                                           [KroneckerDelta(i, j) for
                                                                                            (i, j) in Indices], 1)
            elif Sequence == (1, 0):
                return Rational(1, 2) * -I * (f1[1] + f1[1].conjugate()) * functools.reduce(operator.mul,
                                                                                            [KroneckerDelta(i, j) for
                                                                                             (i, j) in Indices], 1)
            elif Sequence == (1, 1):
                return Rational(1, 2) * (f1[1] - f1[1].conjugate()) * functools.reduce(operator.mul,
                                                                                       [KroneckerDelta(i, j) for (i, j)
                                                                                        in Indices], 1)

    def C2(self, hypercharge):
        """In the Case of U1 the irrep cannot just be F or S, it has to be more precise, for instance it has to match one of the entry of the Particle dictionary defined in the model class."""
        if hypercharge == 'G':
            return 0
        else:
            try:
                return Rational(hypercharge) ** 2
            except TypeError as err:
                return hypercharge**2


    def S2(self, hypercharge):
        """The dynkin index is equal to the Casimir operator for U1 i.e. the Y**2"""
        if hypercharge == 'G':
            return 0
        else:
            try:
                return Rational(hypercharge) ** 2
            except TypeError as err:
                return hypercharge**2

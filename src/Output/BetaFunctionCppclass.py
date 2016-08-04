from ctypes import *


# Wrapper to load c++ betafunction and use them with RGE class in order to solve them, see SolvingCpp.py

class BetaFunctionCpp:
    def __init__(self, namelib, dim):
        self.lib = cdll.LoadLibrary('{}'.format(namelib))
        # dim of the system
        self.dim = dim

    def betafunction(self, t, y, Assumptions={'two-loop': False, 'diag': True}):
        T = c_double(t)
        Y = (c_double * self.dim)()
        F = (c_double * self.dim)()
        TL = c_bool(Assumptions['two-loop'])
        if 'kin2l' in Assumptions:
            Kin2L = c_bool(Assumptions['kin2l'])
        for i in range(self.dim):
            Y[i] = y[i]
        if 'kin2l' in Assumptions:
            self.lib.beta_function(T, byref(Y), byref(F), TL, Kin2L)
        else:
            self.lib.beta_function(T, byref(Y), byref(F), TL)
        return [el for el in F]

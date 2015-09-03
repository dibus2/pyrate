#!/usr/bin/env python 
from RGEsDefinition import *
from ModelsClass import Model, loggingInfo, loggingDebug, loggingCritical
from RGEsmathModule import DynkinIndex, DynkinCasimir, multiplicity, SimplifyTraces, DeterminOrdering, Tr


###########################
# Fermion Anomalous dimension
###########################


def set_globalsFermionAnomalous(model):
    """coppy the globals from RGEsModule into the local gloabls()"""
    for key, val in model.glob.items():
        globals()[key] = val


def CompileFermionAnomalous(Translated, mqt, comb, model, Weyl, GutNorm):
    LH = model.FermionAnomalousToCalculate[mqt][1]
    FinalBeta = (sum([el[0] for el in Translated[mqt]]).doit() / LH).expand()
    if Weyl:
        FinalBeta = FinalBeta.subs(kappa, Rational(1, 2))
    else:
        FinalBeta = FinalBeta.subs(kappa, 1)
    if GutNorm:
        FinalBeta = FinalBeta.subs(model.UGaugeGroups[0][1].g, sqrt(Rational(3, 5)) * model.UGaugeGroups[0][1].g)
    FinalBeta = FinalBeta.subs(tuple([(el, 0) for el in ListAllSymbols['FermionAnomalous']]))
    FinalBeta = DeterminOrdering(model, FinalBeta)
    return FinalBeta


def CFermionAnoI(powe, comb, model):
    """Calculates Eq 28"""
    assert powe[1] == 0
    f1, f2 = comb
    res = Rational(1, 2) * model.Expand(((_Y, s1, f1, p2), (_Ya, s1, p2, f2)))
    powe[0] = powe[0].subs(FermionAnoI, res)
    return powe[0]


def CFermionAnoIg2(powe, comb, model):
    """Calculates Eq 28"""
    assert powe[1] == 2
    f1, f2 = comb
    conserv = functools.reduce(operator.mul, [KroneckerDelta(a, b) for a, b in zip(f1, f2)])
    if conserv != 0:
        res = model.Expand(((_G, gg1), (_C, gg1, f1))) * conserv
    else:
        res = 0
    powe[0] = powe[0].subs(FermionAnoIg2, res)
    return powe[0]

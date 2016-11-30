#!/usr/bin/env python 
from RGEsDefinition import *
from ModelsClass import Model, loggingInfo, loggingDebug, loggingCritical
from RGEsmathModule import DynkinIndex, DynkinCasimir, multiplicity, SimplifyTraces, DeterminOrdering, Tr


###########################
# Scalar Anomalous dimension
###########################


def set_globalsScalarAnomalous(model):
    """coppy the globals from RGEsModule into the local gloabls()"""
    for key, val in model.glob.items():
        globals()[key] = val


def CompileScalarAnomalous(Translated, mqt, comb, model, Weyl, GutNorm):
    LH = model.ScalarAnomalousToCalculate[mqt][1]
    FinalBeta = (sum([el[0] for el in Translated[mqt]]).doit() / LH).expand()
    if Weyl:
        FinalBeta = FinalBeta.subs(kappa, Rational(1, 2))
    else:
        FinalBeta = FinalBeta.subs(kappa, 1)
    if GutNorm:  # Apply any GUT normalization
        FinalBeta = FinalBeta.subs(model.GutNorm)
    FinalBeta = FinalBeta.subs(tuple([(el, 0) for el in ListAllSymbols['ScalarAnomalous']]))
    FinalBeta = DeterminOrdering(model, FinalBeta)
    return FinalBeta


def CY2ab(powe, comb, model):
    """Calculate the Yab term."""
    assert powe[1] == 0
    sc1, sc2 = comb
    res = model.Expand((('Yab2S', sc1, sc2)))
    if res != 0:
        res = res.doit()
    powe[0] = powe[0].subs(Y2ab, res)
    return powe[0]


def CScalarAnoI(powe, comb, model):
    """Calculate the ScalarAnoI term."""
    assert powe[1] == 2
    sc1, sc2 = comb
    conserv = functools.reduce(operator.mul, [KroneckerDelta(a, b) for a, b in zip(sc1, sc2)])
    if conserv == 1:
        res = (3 - xi) * model.Expand(((_G, gg1), (_Cs, gg1, sc1)))
    else:
        res = 0
    powe[0] = powe[0].subs(ScalarAnoI, res)
    return powe[0]


########################
# Two-loop contributions
########################

def CL2ab(powe, comb, model):
    """Calculates the term L2ab Eq 24"""
    assert powe[1] == 0
    sc1, sc2 = comb
    res = Rational(1, 6) * model.Expand(((_L, sc1, s1, s2, s3), (_L, sc2, s1, s2, s3)))
    powe[0] = powe[0].subs(L2ab, res)
    return powe[0]


def CH2ab(powe, comb, model):
    """Calculates the term Hab Eq 25"""
    assert powe[1] == 0
    sc1, sc2 = comb
    res = Rational(1, 2) * (model.Expand((('Chain3Y', sc1, sc2, s1, p1, p2), (_Ya, s1, p2, p1)))
                            + model.Expand((('Chain3Ya', sc1, sc2, s1, p1, p2), (_Y, s1, p2, p1)))
                            )
    powe[0] = powe[0].subs(H2ab, res)
    return powe[0]


def CHbar2ab(powe, comb, model):
    """Calculates the term Hbarab Eq 26"""
    assert powe[1] == 0
    sc1, sc2 = comb
    res = Rational(1, 2) * (model.Expand((('Chain3Y', sc1, s1, sc2, p1, p2), (_Ya, s1, p2, p1)))
                            + model.Expand((('Chain3Ya', sc1, s1, sc2, p1, p2), (_Y, s1, p2, p1)))
                            )
    powe[0] = powe[0].subs(Hbar2ab, res)
    return powe[0]


def CY2Fab(powe, comb, model):
    """Calculates the term CY2Fab Eq 27"""
    assert powe[1] == 2
    sc1, sc2 = comb
    res = Rational(1, 2) * (model.Expand2(
        ((_G, gg1), (_C, gg1, p1), [[(_Y, sc1, p1, p2), (_Ya, sc2, p2, p1)], [(_Y, sc2, p1, p2), (_Ya, sc1, p2, p1)]]))
                            )
    powe[0] = powe[0].subs(Y2Fab, res)
    return powe[0]


def CC2SC2S(powe, comb, model):
    """Calculates the term C_2(S)^ac C_2(S) ^cb Eq 23"""
    assert powe[1] == 4
    sc1, sc2 = comb
    conserv = functools.reduce(operator.mul, [KroneckerDelta(a, b) for a, b in zip(sc1, sc2)])
    if conserv == 1:
        res = model.Expand(((_G, gg1), (_G, gg2), (_Cs, gg1, sc1), (_Cs, gg2, sc1)))
    else:
        res = 0
    powe[0] = powe[0].subs(C2SC2S, res)
    return powe[0]


def CScalarAnoII(powe, comb, model):
    """Calculates the term C_2(S)^ab (35/3-2xi ...) Eq 23"""
    assert powe[1] == 4
    sc1, sc2 = comb
    conserv = functools.reduce(operator.mul, [KroneckerDelta(a, b) for a, b in zip(sc1, sc2)])
    if conserv == 1:
        res = ((Rational(35, 3) - 2 * xi - Rational(1, 4) * xi ** 2) * model.Expand(
            ((_G, gg1), (_G, gg1), (_Cs, gg1, sc1), (_Cg, gg1)))
               - Rational(10, 3) * kappa * model.Expand(((_G, gg1), (_G, gg1), (_Cs, gg1, sc1), (_S, gg1)))
               - Rational(11, 12) * model.Expand(((_G, gg1), (_G, gg1), (_Cs, gg1, sc1), (_Ss, gg1)))
               )
    else:
        res = 0
    powe[0] = powe[0].subs(ScalarAnoII, res)
    return powe[0]

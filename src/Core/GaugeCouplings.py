#!/usr/bin/env python 
import functools, operator, copy
import re as regular  # avoid the double definition of real part of sympy
from RGEsDefinition import *
from ModelsClass import Model, loggingInfo, loggingDebug, loggingCritical
from RGEsmathModule import DynkinIndex, DynkinCasimir, DeterminOrdering, summatrix, multiplicityKin, multiplicity


################################################################################
# Gauge Couplings transformation rules including kinetic mixing
################################################################################

# powe 3
def CompileGaugeCouplings(model, Translated, x, Weyl, GutNorm):
    if not (x == 'abelian'):
        FinalBeta = sum([el[0] for el in Translated[x]]).expand()
    else:
        FinalBeta = summatrix([el[0] for el in Translated[x]]).expand()
    if Weyl:
        FinalBeta = FinalBeta.subs(kappa, Rational(1, 2))
    else:
        FinalBeta = FinalBeta.subs(kappa, 1)
    if GutNorm:  # Apply any GUT normalization
        FinalBeta = FinalBeta.subs(model.GutNorm)
        if model.kinmixing and x == 'abelian':
            FinalBeta = Matrix([el * model.facGutNorm[iel] for iel, el in enumerate(FinalBeta)]).reshape(2, 2).expand()
        else:
            if model.GetGroupFromName[x].U:
                # find in which position is the corresponding factor
                fac = [model.facGutNorm[iel] for iel,el in enumerate(model.UsectorMatrix) if el == model.GetGroupFromName[x].g]
                if fac and len(fac) == 1:
                    FinalBeta = (FinalBeta * fac[0]).expand()
                elif fac and len(fac) != 1:
                    loggingCritical("Error in applying Gut normalization, contact the authors", verbose=True)
    FinalBeta = DeterminOrdering(model, FinalBeta)
    return FinalBeta


def set_globals(model):
    """coppy the globals from RGEsModule into the local gloabls()"""
    for key, val in model.glob.items():
        globals()[key] = val


def abeliansector1loop(model):
    """returns the 1-Loop contribution for the abelian sector when kinetic mixing is in"""
    out = [Rational(4, 3) * kappa * model.UsectorMatrix * summatrix(
        [pp.W * pp.W.transpose() * multiplicityKin(pp, model) for pp in model.Fermions.values()]),
           Rational(1, 6) * model.UsectorMatrix * summatrix(
               [pp.W * pp.W.transpose() * multiplicityKin(pp, model) for pp in model.Scalars.values()])]
    return 1 / (4 * pi) ** 2 * (out[0] + out[1]).expand()


def ApplyEq106(powe, Mod, x):
    """Apply the transformation rule defined by Eq. 106
	No modification for the kinetic mixing which means that they must be 
	multiplied directly by G^3"""
    assert powe[1] == 3
    if not (x[1].U):
        return powe[0].subs([(C2(G), x[1].C2('G'))]) if powe[1] == 3 else powe[0]
    else:
        "There is no Casimir for the adjoint for the U1 group"
        return powe[0].subs([(C2(G), 0)]) if powe[1] == 3 else powe[0]


def ApplyEq107(powe, Model, x):
    """Apply the transformation rule defined by Eq. 107"""
    assert powe[1] == 3
    out = powe[0].subs([(S2(F), DynkinIndex(Model, x, 'Fermions')), (S2(S), DynkinIndex(Model, x, 'Scalars'))])
    return out


# power 5
def abeliansector2loop(model):
    """Returns the g^5 contribution for the abelian sector, Eq.(14) of 1308.1674"""
    out = [model.UsectorMatrix * summatrix([pp.W * pp.W.transpose() * multiplicityKin(pp, model) * (
        sum([g.g ** 2 * g.C2(pp.Qnb[name]) for name, g, t in model.NonUGaugeGroups]) + (pp.W.transpose() * pp.W)[0])
                                            for pp in Vals])
           for Vals in [model.Fermions.values(), model.Scalars.values()]]
    # coefficient from Eq.(30)
    if out[0] != 0:
        res = 1 / (4 * pi) ** 4 * (4 * kappa * out[0])
    if out[1] != 0:
        res += 1 / (4 * pi) ** 4 * (2 * out[1])
    outY4F = -1 / (4 * pi) ** 4 * 2 * kappa * model.UsectorMatrix * model.Expand(
        (('WW_T', p1), (_Y, s1, p1, p2), (_Ya, s1, p2, p1)))
    if outY4F != 0:
        res += outY4F
    return res


def ApplyEq108(powe, Mod, x):
    """Apply the transformation rule defined by Eq. 108."""
    assert powe[1] == 5
    if not (x[1].U):
        loggingDebug('Applying equation 108 on : {}'.format(powe))
        return powe[0].subs(C2(G) ** 2, (x[1].C2('G')) ** 2)
    else:
        return powe[0].subs(C2(G), 0)


def ApplyEq109(powe, Mod, x):
    """Apply the transformation rule defined by Eq. 109."""
    assert powe[1] == 5
    if not (x[1].U):
        return powe[0].subs([(C2(G) * S2(F), x[1].C2('G') * DynkinIndex(Mod, x, 'Fermions')),
                             (C2(G) * S2(S), x[1].C2('G') * DynkinIndex(Mod, x, 'Scalars'))])
        loggingDebug('Applying equation 109 on : {}'.format(powe))
    else:
        return powe[0].subs([(C2(G) * S2(F), 0), (C2(G) * S2(S), 0)])


def ApplyEq110(powe, model, x):
    """Apply the transformation rule defined by Eq. 110."""
    assert powe[1] == 5
    if model.kinmixing:
        # if kinetic mixing apply Eq. 16 of FS
        powe[0] = powe[0].expand().subs([(x[1].g ** 5 * C2(F) * S2(F), x[1].g ** 3 * sum([
                                                                                             x[1].S2(pp.Qnb[x[
                                                                                                 0]]) * multiplicity(x,
                                                                                                                     '',
                                                                                                                     pp,
                                                                                                                     model) * (
                                                                                             sum(
                                                                                                 [g.g ** 2 * g.C2(
                                                                                                     pp.Qnb[name]) for
                                                                                                  name, g, t in
                                                                                                  model.NonUGaugeGroups]) +
                                                                                             (pp.W.transpose() * pp.W)[
                                                                                                 0]) for pp in
                                                                                             model.Fermions.values()])),
                                         # and the scalars identical
                                         (x[1].g ** 5 * C2(S) * S2(S), x[1].g ** 3 * sum([
                                                                                             x[1].S2(pp.Qnb[x[
                                                                                                 0]]) * multiplicity(x,
                                                                                                                     '',
                                                                                                                     pp,
                                                                                                                     model) * (
                                                                                             sum(
                                                                                                 [g.g ** 2 * g.C2(
                                                                                                     pp.Qnb[name]) for
                                                                                                  name, g, t in
                                                                                                  model.NonUGaugeGroups]) +
                                                                                             (pp.W.transpose() * pp.W)[
                                                                                                 0]) for pp in
                                                                                             model.Scalars.values()]))])
        return powe[0]
    else:
        powe[0] = powe[0].expand().subs([(x[1].g ** 5 * C2(F) * S2(F), DynkinCasimir(model, x, 'Fermions')),
                                         (x[1].g ** 5 * C2(S) * S2(S), DynkinCasimir(model, x, 'Scalars'))])
    return powe[0]


def Y4F(powe, model, x):
    """Calculates the yukawa term in the gauge coupling, i.e. Y4(F)."""
    assert powe[1] == 3
    # Expand the product 1/d*C2(i)*Y(z,i,j)*Ydagger(z,j,i) z among higgs fields and i,j among fermions
    res = Rational(1, x[1].dimAdj) * model.Expand(((_C, x[0], p1), (_Y, s1, p1, p2), (_Ya, s1, p2, p1)))
    return powe[0].expand().subs(Y4(F), res)

#!/usr/bin/env python 
from RGEsDefinition import *
from RGEsmathModule import SimplifyTraces, DeterminOrdering


##################################################################################
# Scalar mass functions
##################################################################################


def set_globalsScM(model):
    """coppy the globals from RGEsModule into the local gloabls()"""
    for key, val in model.glob.items():
        globals()[key] = val


def CompileScalar(Translated, mqt, comb, model, Weyl, GutNorm):
    LH = model.ScMToCalculate[mqt][1]
    FinalBeta = (sum([el[0] for el in Translated[mqt]]).doit() / LH).expand()
    if Weyl:
        FinalBeta = FinalBeta.subs(kappa, Rational(1, 2))
    else:
        FinalBeta = FinalBeta.subs(kappa, 1)
    if GutNorm:
        FinalBeta = FinalBeta.subs(model.UGaugeGroups[0][1].g, sqrt(Rational(3, 5)) * model.UGaugeGroups[0][1].g)
    FinalBeta = FinalBeta.subs(tuple([(el, 0) for el in ListAllSymbols['ScalarMasses']]))
    FinalBeta = DeterminOrdering(model, FinalBeta)
    return FinalBeta


######
# 1Loop
######

def Cmeflabef(powe, comb, model):
    """Calculates the first term of Eq (86)"""
    sc1, sc2 = comb
    res = model.Expand(((_ms, s1, s2), (_L, sc1, sc2, s1, s2)))
    powe[0] = powe[0].subs(meflabef, res)
    return powe[0]


def Chaefhbef(powe, comb, model):
    """Calculates the second term in Eq(86)"""
    if model.TriToCalculate != {}:
        sc1, sc2 = comb
        res = model.Expand(((_h, sc1, s1, s2), (_h, sc2, s1, s2)))
        powe[0] = powe[0].subs(haefhbef, res)
    else:
        powe[0] = powe[0].subs(haefhbef, 0)
    return powe[0]


def CLSab(powe, comb, model):
    """Calculate the term LSab Eq 88"""
    sc1, sc2 = comb
    res = [model.Expand(((_G, gg1), (_Cs, gg1, i), (_ms, sc1, sc2)))
           for i in [sc1, sc2]
           ]
    res = sum(res)
    reskin = 0
    if model.kinmixing:
        reskin = sum([model.Expand2(((_Ckins,i), (_ms, sc1, sc2))) for i in [sc1, sc2]])
        if reskin != 0:
            reskin = reskin.doit()
    if res != 0:
        res = res.doit()
    res += reskin
    powe[0] = powe[0].subs(LSab, res)
    return powe[0]


def CLYab(powe, comb, model):
    """Calculate the LYabcd term Eq 89"""
    sc1, sc2 = comb
    res = [model.Expand((('Yab2S', sc, sc), (_ms, sc1, sc2)))
           for sc in comb]
    res = sum(res)
    if res != 0:
        res = res.doit()
    powe[0] = powe[0].subs(LYab, res)
    return powe[0]


def CHab(powe, comb, model):
    """Calculates the Habcd term Eq 87"""
    if model.FMToCalculate != {}:
        sc1, sc2 = comb
        res = SimplifyTraces(model.Expand((('Hab', sc1, sc2))), model)
        powe[0] = powe[0].subs(Hab, res)
    else:
        powe[0] = powe[0].subs(Hab, 0)
    return powe[0]


#######
# 2-Loop
#######

def CL2mab(powe, comb, model):
    """Calculates the term L2mab Eq 90 line 1 first term"""
    sc1, sc2 = comb
    res = [Rational(1, 2) * model.Expand((('L2abS', a, a), (_ms, sc1, sc2)), dotrace=False)
           for a in [sc1, sc2]
           ]
    res = sum(res)
    if res != 0:
        res = res.doit()
    powe[0] = powe[0].subs(L2mab, res)
    return powe[0]


def CLbar3ab(powe, comb, model):
    """Calculates the term Lbar3ab Eq 91"""
    a, b = comb
    res = (model.Expand(((_L, a, b, s1, s2), (_h, s1, s3, s4), (_h, s2, s3, s4)))
           + 2 * model.Expand(((_ms, s1, s2), (_L, a, s1, s3, s4), (_L, b, s2, s3, s4)))
           + 2 * model.Expand(((_h, a, s1, s2), (_h, s2, s3, s4), (_L, b, s1, s3, s4)))
           + 2 * model.Expand(((_h, b, s1, s2), (_h, s2, s3, s4), (_L, a, s1, s3, s4)))
           )
    powe[0] = powe[0].subs(Lbar3ab, res)
    return powe[0]


def CLbar2Yab(powe, comb, model):
    """Calculates the term Lbar2Yab Eq 92"""
    a, b = comb
    res = (model.Expand((('Yab2S', s1, s2), (_ms, s3, s2), (_L, a, b, s3, s1)))
           + model.Expand((('Yab2S', s1, s2), (_h, a, s3, s1), (_h, b, s3, s2)))
           )
    powe[0] = powe[0].subs(Lbar2Yab, res)
    return powe[0]


def CHbarlab(powe, comb, model):
    """Calculates the term Hbarlab Eq 93"""
    a, b = comb
    res1 = Rational(1, 2) * (
        model.Expand2(((_L, a, b, s1, s2), ([[(_mf, p1, p2), (_Ya, s1, p2, p3), (_mf, p3, p4), (_Ya, s2, p4, p1)],
                                             [(_Y, s2, p1, p2), (_mfa, p2, p3), (_Y, s1, p3, p4), (_mfa, p4, p1)]])))
        # line 1
    )
    res2 = model.Expand2(((_ms, s1, s2), ([[('Chain3Y', a, s1, b, p1, p2), (_Ya, s2, p2, p1)],
                                           [('Chain3Y', s2, b, s1, p1, p2), (_Ya, a, p2, p1)]])))  # Line 2

    res3 = model.Expand2(((_h, a, s1, s2), ([[('Chain3Ya', s2, b, s1, p1, p2), (_mf, p2, p1)],
                                             [(_mfa, p1, p2), ('Chain3Y', s1, b, s2, p2, p1)]])))  # Line 3

    res4 = model.Expand2(((_h, b, s1, s2), ([[('Chain3Ya', s2, a, s1, p1, p2), (_mf, p2, p1)],
                                             [(_mfa, p1, p2), ('Chain3Y', s1, a, s2, p2, p1)]])))  # Line 4
    resall = res1 + res2 + res3 + res4
    powe[0] = powe[0].subs(Hbarlab, resall)
    return powe[0]


def CH2mab(powe, comb, model):
    """Calculates H2(i)mab term Eq 90 line 2"""
    a, b = comb
    res = [model.Expand((('H2abS', i, i), (_ms, a, b)))
           for i in comb
           ]
    res = sum(res)
    powe[0] = powe[0].subs(H2mab, res)
    return powe[0]


def CHbar2mab(powe, comb, model):
    """Calculates H2bar2mab term Eq 90 line 2"""
    a, b = comb
    res = [model.Expand((('Hbar2abS', i, i), (_ms, a, b)), dotrace=False)
           for i in comb
           ]
    res = sum(res)
    powe[0] = powe[0].subs(Hbar2mab, res)
    return powe[0]


def CHYabHbarYabH3ab(powe, comb, model):
    """Calculates HYab and HYbarab Eq 94 and 93 and 96."""
    a, b = comb
    resHY = (
        2 * model.Expand2(
            (([[('Y2F', p1, p2), (_mfa, p2, p3), (_mf, p3, p4)], [(_mfa, p1, p2), (_mf, p2, p3), ('Y2F', p3, p4)]]),
             ([[(_Ya, a, p4, p5), (_Y, b, p5, p1)], [(_Ya, b, p4, p5), (_Y, a, p5, p1)]])))  # Line 1 HY

        + 2 * model.Expand2(
            (('Y2F', p1, p2), ([[(_Ya, a, p2, p3), (_mf, p3, p4)], [(_mfa, p2, p3), (_Y, a, p3, p4)]]),
             ([[(_Ya, b, p4, p5), (_mf, p5, p1)], [(_mfa, p4, p5), (_Y, b, p5, p1)]])))  # Line 2 and Line 3
        + 2 * model.Expand2(
            (('Y2F', p1, p2), ([[(_Ya, b, p2, p3), (_mf, p3, p4)], [(_mfa, p2, p3), (_Y, b, p3, p4)]]),
             ([[(_Ya, a, p4, p5), (_mf, p5, p1)], [(_mfa, p4, p5), (_Y, a, p5, p1)]])))  # Line 4 and 5
    )
    resHYbar = 2 * (model.Expand2(
        (([[(_Y, s1, p1, p2), ('Chain3Ya', a, s1, b, p2, p3)], [(_Y, s1, p1, p2), ('Chain3Ya', b, s1, a, p2, p3)]]),
         (_mf, p3, p4), (_mfa, p4, p1)))  # Line 1

                    + model.Expand2(((_mfa, p1, p2), (_Y, s1, p2, p3), (_mfa, p3, p4),
                                     ([[('Chain3Y', a, b, s1, p4, p1)], [('Chain3Y', b, a, s1, p4, p1)]])))  # Line 2

                    + model.Expand2(((_Y, s1, p1, p2), (
        [[(_Ya, a, p2, p3), (_Y, s1, p3, p4), (_mfa, p4, p5)], [(_mfa, p2, p3), (_Y, s1, p3, p4), (_Ya, a, p4, p5)]]),
                                     ([[(_Y, b, p5, p6), (_mfa, p6, p1)],
                                       [(_mf, p5, p6), (_Ya, b, p6, p1)]])))  # Line 3 and 4

                    + model.Expand2(((_Y, s1, p1, p2), (
        [[(_Ya, b, p2, p3), (_Y, s1, p3, p4), (_mfa, p4, p5)], [(_mfa, p2, p3), (_Y, s1, p3, p4), (_Ya, b, p4, p5)]]),
                                     ([[(_Y, a, p5, p6), (_mfa, p6, p1)], [(_mf, p5, p6), (_Ya, a, p6, p1)]])))
                    # Line 5 and 6
                    )
    resH3ab = 2 * (model.Expand2((([[('Chain3Y', a, b, s1, p1, p2)], [('Chain3Y', b, a, s1, p1, p2)]]),
                                  (_mfa, p2, p3), (_mf, p3, p4), (_Ya, s1, p4, p1)))  # Line 1

                   + model.Expand2(((_Ya, s1, p1, p2), (_mf, p2, p3), (_mfa, p3, p4),
                                    ([[('Chain3Y', s1, a, b, p4, p1)], [('Chain3Y', s1, b, a, p4, p1)]])))  # Line 2

                   + model.Expand2((([[(_Y, a, p1, p2), (_mfa, p2, p3)], [(_mf, p1, p2), (_Ya, a, p2, p3)]]),
                                    ([[(_Y, s1, p3, p4), (_Ya, b, p4, p5), (_mf, p5, p6), (_Ya, s1, p6, p1)],
                                      [(_Y, s1, p3, p4),
                                       (_mfa, p4, p5), (_Y, b, p5, p6), (_Ya, s1, p6, p1)]])))  # Line 3 and 4

                   + model.Expand2((([[(_Y, b, p1, p2), (_mfa, p2, p3)], [(_mf, p1, p2), (_Ya, b, p2, p3)]]),
                                    ([[(_Y, s1, p3, p4), (_Ya, a, p4, p5), (_mf, p5, p6), (_Ya, s1, p6, p1)],
                                      [(_Y, s1, p3, p4), (_mfa, p4, p5), (_Y, a, p5, p6), (_Ya, s1, p6, p1)]])))
                   # Line 5
                   )
    resHY = resHY.doit() if resHY != 0 else resHY
    resHYbar = resHYbar.doit() if resHYbar != 0 else resHYbar
    resH3ab = resH3ab.doit() if resH3ab != 0 else resH3ab
    resHYbarad = resHYbar.adjoint() if resHYbar != 0 else 0
    resHYbarad = resHYbarad.doit() if resHYbarad != 0 else 0
    resFinal = SimplifyTraces(resHY + resHYbar + resHYbarad + resH3ab, model)
    powe[0] = powe[0].subs(HYab, resFinal).subs(HbarYab, 0).subs(H3ab, 0)
    return powe[0]


#####
# g**2
#####

def CLbar2Sab(powe, comb, model):
    """Calculates the term Lbar2Sab Eq 97"""
    a, b = comb
    res = (model.Expand(((_G, gg1), (_Cs, gg1, s1), (_L, a, b, s2, s1), (_ms, s2, s1)))
           + model.Expand(((_G, gg1), (_Cs, gg1, s1), (_h, a, s2, s1), (_h, b, s2, s1)))
           )
    reskin = 0
    if model.kinmixing:
        reskin = (model.Expand(((_Ckins, s1), (_L, a, b, s2, s1), (_ms, s2, s1)))
        + model.Expand(((_Ckins, s1), (_h, a, s2, s1), (_h, b, s2, s1)))
                  )
    res += reskin
    powe[0] = powe[0].subs(Lbar2Sab, res)
    return powe[0]


def CL2gab(powe, comb, model):
    """Calculates term L2gab Eq 98"""
    a, b = comb
    res = (model.Expand(((_G, gg1), (_L, a, b, s1, s2), (_ms, s3, s4), (_Th, gg1, s1, s3), (_Th, gg1, s2, s4)),
                        MatStruc=['a', 'a'])
           + model.Expand(((_G, gg1), (_h, a, s1, s2), (_h, b, s3, s4), (_Th, gg1, s1, s3), (_Th, gg1, s2, s4)),
                          MatStruc=['a', 'a'])
           )
    reskin = 0
    if model.kinmixing:
        reskin = (
            model.Expand(((_L, a, b, s1, s2), (_ms, s3, s4), ('Thetakin', s1, s3, s2, s4)))
            + model.Expand(((_h, a, s1, s2), (_h, b, s3, s4), ('Thetakin', s1, s3, s2, s4)))
            )
    res += reskin
    powe[0] = powe[0].subs(L2gab, res)
    return powe[0]


def CY2Fmab(powe, comb, model):
    """Calculates the Y2Fab term"""
    a, b = comb
    res = [model.Expand(((_G, gg1), ('Y2FabS', gg1, i, i), (_ms, a, b)))
           for i in comb
           ]
    res = sum(res)
    reskin = 0
    if model.kinmixing:
        reskin = [model.Expand((('Y2FabSkin', i, i), (_ms, a, b)))
                  for i in comb
                  ]
        reskin = sum(reskin)
    res += reskin
    powe[0] = powe[0].subs(Y2Fmab, res)
    return powe[0]


def CHSab(powe, comb, model):
    """Calculates the term CHSab Eq 99"""
    a, b = comb
    res = [model.Expand(((_G, gg1), (_Cs, gg1, i), ('Hab', a, b)))
           for i in comb
           ]
    res = sum(res)
    reskin = 0
    if model.kinmixing:
        reskin = [model.Expand(((_Ckins, i), ('Hab', a, b))) for i in comb]
        reskin = sum(reskin)
    res+=reskin
    powe[0] = powe[0].subs(HSab, res)
    return powe[0]


def CHFab(powe, comb, model):
    """Calculates the term HFab Eq 100"""
    a, b = comb
    res = 2 * (
        model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_Y, a, p1, p2)], [(_Y, a, p1, p2), (_C, gg1, p2)]]),
                       (_Ya, b, p2, p3), (_mf, p3, p4), (_mfa, p4, p1)))  # Line 1

        + model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_Y, b, p1, p2)], [(_Y, b, p1, p2), (_C, gg1, p2)]]),
                         (_Ya, a, p2, p3), (_mf, p3, p4), (_mfa, p4, p1)))  # Line 2

        + model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_C, gg1, p2)]]),
                         (_mfa, p2, p3),
                         ([[(_Y, a, p3, p4), (_Ya, b, p4, p1)], [(_Y, b, p3, p4), (_Ya, a, p4, p1)]])))  # Line 3

        + model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_Y, a, p1, p2)], [(_Y, a, p1, p2), (_C, gg1, p2)]]),
                         (_mfa, p2, p3),
                         ([[(_Y, b, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, b, p4, p1)]])))  # Line 4

        + model.Expand2(
            ((_G, gg1), ([[(_C, gg1, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_C, gg1, p2)]]), (_Ya, a, p2, p3),
             ([[(_Y, b, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, b, p4, p1)]])))  # Line 5

        + model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_Y, b, p1, p2)], [(_Y, b, p1, p2), (_C, gg1, p2)]]),
                         (_mfa, p2, p3),
                         ([[(_Y, a, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, a, p4, p1)]])))  # Line 6

        + model.Expand2(((_G, gg1), ([[(_C, gg1, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_C, gg1, p2)]]),
                         (_Ya, b, p2, p3), ([[(_Y, a, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, a, p4, p1)]])))
        # Line 7
    )
    reskin = 0
    if model.kinmixing:
        reskin = 2 * (
              model.Expand2((([[(_Ckin, p1), (_Y, a, p1, p2)], [(_Y, a, p1, p2), (_Ckin, p2)]]),
                                     (_Ya, b, p2, p3), (_mf, p3, p4), (_mfa, p4, p1)))
            + model.Expand2((([[(_Ckin, p1), (_Y, b, p1, p2)], [(_Y, b, p1, p2), (_Ckin, p2)]]),
                             (_Ya, a, p2, p3), (_mf, p3, p4), (_mfa, p4, p1)))
            + model.Expand2((([[(_Ckin, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_Ckin, p2)]]),
                         (_mfa, p2, p3),([[(_Y, a, p3, p4), (_Ya, b, p4, p1)], [(_Y, b, p3, p4), (_Ya, a, p4, p1)]])))
            + model.Expand2((([[(_Ckin, p1), (_Y, a, p1, p2)], [(_Y, a, p1, p2), (_Ckin, p2)]]),
                         (_mfa, p2, p3),
                         ([[(_Y, b, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, b, p4, p1)]])))
            + model.Expand2((([[(_Ckin, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_Ckin, p2)]]), (_Ya, a, p2, p3),
             ([[(_Y, b, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, b, p4, p1)]])))
            + model.Expand2((([[(_Ckin, p1), (_Y, b, p1, p2)], [(_Y, b, p1, p2), (_Ckin, p2)]]),
                         (_mfa, p2, p3),
                         ([[(_Y, a, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, a, p4, p1)]])))  # Line 6
            + model.Expand2((([[(_Ckin, p1), (_mf, p1, p2)], [(_mf, p1, p2), (_Ckin, p2)]]),
                         (_Ya, b, p2, p3), ([[(_Y, a, p3, p4), (_mfa, p4, p1)], [(_mf, p3, p4), (_Ya, a, p4, p1)]])))
        )
    res += reskin
    powe[0] = powe[0].subs(HFab, res)
    return powe[0]


#####
# g**4
#####

def CLSabg4(powe, comb, model):
    """Calculates the term LSabg4 Eq 90 3rd line from the bottom"""
    a, b = comb
    res = [Rational(35, 3) * model.Expand(((_G, gg1), (_G, gg1), (_Cg, gg1), (_Cs, gg1, i), (_ms, a, b)))
           - Rational(10, 3) * kappa * model.Expand(((_G, gg1), (_G, gg1), (_S, gg1), (_Cs, gg1, i), (_ms, a, b)))
           - Rational(11, 12) * model.Expand(((_G, gg1), (_G, gg1), (_Ss, gg1), (_Cs, gg1, i), (_ms, a, b)))
           for i in comb
           ]
    res = sum(res)
    reskin = 0
    if model.kinmixing:
        reskin = [
            - Rational(10, 3) * kappa * model.Expand(((_SfCkins, i), (_ms, a, b)))
            - Rational(11, 12) * model.Expand(((_SsCkins, i), (_ms, a, b)))
            for i in comb
        ]
        reskin = sum(reskin)
    res += reskin
    powe[0] = powe[0].subs(LSabg4, res)
    return powe[0]


def CLSSab(powe, comb, model):
    """Calculates the term LSabg4 Eq 90 before last line"""
    a, b = comb
    res = [model.Expand(((_G, gg1), (_G, gg2), (_Cs, gg1, i), (_Cs, gg2, i), (_ms, a, b)))
           for i in comb
           ]
    res = sum(res)
    reskin = 0
    if model.kinmixing:
        reskin = [
            2 * model.Expand(((_G,gg1), (_Cs, gg1, i), (_Ckins, i), (_ms, a, b)))# 2* C Ckin
            + model.Expand(((_Ckins, i), (_Ckins, i), (_ms, a, b)))#Ckin Ckin
        ]
        reskin = sum(reskin)
    res += reskin
    powe[0] = powe[0].subs(LSSab, res)
    return powe[0]


def CAlab(powe, comb, model):
    """Calculates the term Alab Eq 90 before last line """
    a, b = comb
    res = 2 * model.Expand(((_G, gg1), (_G, gg2), (_ms, s1, s2), ('Theta4', gg1, gg2, s1, s2, a, b)))
    reskin = 0
    if model.kinmixing:
        reskin = 2 * (
            model.Expand(((_ms, s1, s2), ('Theta2g2kin', s1, a, s3, s4),
                      ('Thetakin', s3, s2, s4, b)))
        + model.Expand(((_ms, s1, s2), ('Theta2g2kin', s3, s4, s2, b),
                        ('Thetakin', s1, s3, a, s4)))
        + model.Expand(((_ms, s1, s2), ('Theta2g2kin', s1, s4, s3, b),
                        ('Thetakin', s3, s2, a, s4)))
        + model.Expand(((_ms, s1, s2), ('Theta2g2kin', s3, a, s2, s4),
                        ('Thetakin', s1, s3, s4, b)))
        # Terms C*C
        + model.Expand(((_ms, s1, s2), ('Thetakin', s1, s3, a, s4),
                        ('Thetakin', s3, s2, s4, b)))
        + model.Expand(((_ms, s1, s2), ('Thetakin', s1, s3, s4, b),
                        ('Thetakin', s3, s2, a, s4)))
        )
    res += reskin
    powe[0] = powe[0].subs(Alab, res)
    return powe[0]


def CAbarlab(powe, comb, model):
    """Calculates the term Alab Eq 90 before last line """
    a, b = comb
    res = 2 * model.Expand(((_G, gg1), (_G, gg2), (_ms, s1, s2), ('Theta4', gg1, gg2, a, s1, b, s2)))
    reskin = 0
    if model.kinmixing:
        reskin = 2 * (
            model.Expand(((_ms, s1, s2), ('Theta2g2kin', a, b, s3, s4), ('Thetakin', s3, s1, s4, s2)))
            + model.Expand(((_ms, s1, s2), ('Theta2g2kin', s3, s4, s1, s2), ('Thetakin', a, s3, b, s4)))
            + model.Expand(((_ms, s1, s2), ('Theta2g2kin', a, s4, s3, s2), ('Thetakin', s3, s1, b, s4)))
            + model.Expand(((_ms, s1, s2), ('Theta2g2kin', s3, b, s1, s4), ('Thetakin', a, s3, s4, s2)))
            # Terms C*C
            + model.Expand(((_ms, s1, s2), ('Thetakin', a, s3, b, s4), ('Thetakin', s3, s1, s4, s2)))
            + model.Expand(((_ms, s1, s2), ('Thetakin', a, s3, s4, s2), ('Thetakin', s3, s1, b, s4)))
        )
    res += reskin
    powe[0] = powe[0].subs(Abarlab, res)
    return powe[0]


def CBYabBbarYab(powe, comb, model):
    """Calculates the two terms BYab and BbarYab Eq 104 and 105"""
    a, b = comb
    res = (
            model.Expand((
            (_G, gg1), (_G, gg2), (_Th, gg1, a, s1), (_Th, gg2, s1, b), (_Ta, gg3, p1, p2), (_Ta, gg4, p2, p3),
            (_mf, p3, p4), (_mfa, p4, p1)), MatStruc=[['a', 'b', 'a', 'b'], ['b', 'a', 'a', 'b']])
        +   model.Expand(((_G, gg1), (_G, gg2), (_Th, gg1, a, s1), (_Th, gg2, s1, b), (_mf, p1, p2), (_T, gg3, p2, p3),
                        (_T, gg4, p3, p4), (_mfa, p4, p1)), MatStruc=[['a', 'b', 'a', 'b'], ['b', 'a', 'a', 'b']])
        -10 * (
            model.Expand(((_G, gg1), (_G, gg2), (_Th, gg1, a, s1), (_Th, gg2, s1, b), (_Ta, gg3, p1, p2), (_mf, p2, p3),
                      (_Ta, gg4, p3, p4), (_mfa, p4, p1)), MatStruc=[['a', 'b', 'a', 'b'], ['b', 'a', 'a', 'b']])
        )
    )
    reskin = 0
    if model.kinmixing:
        reskin = (
            2 * model.Expand((('Thetakin4W', a, p1, a, p1), (_mf, p1, p2), (_mfa, p2, p1)))
            + model.Expand(((_G, gg1),('ThetakinWsWf', a, p1), (_Th, gg1, s1, b), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1),('ThetakinWsWf', a, p2), (_Th, gg1, s1, b), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1),('ThetakinWsWf', b, p1), (_Th, gg1, a, s1), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1),('ThetakinWsWf', b, p2), (_Th, gg1, a, s1), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
            #Y T^A T^B
            + 2 * model.Expand((('Thetakin4W', a, p1, a, p1), (_mf, p2, p1), (_mfa, p1, p2)))
            + model.Expand(((_G, gg1), ('ThetakinWsWf', a, p2), (_Th, gg1, s1, b), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1), ('ThetakinWsWf', a, p3), (_Th, gg1, s1, b), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1), ('ThetakinWsWf', b, p2), (_Th, gg1, a, s1), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
            + model.Expand(((_G, gg1), ('ThetakinWsWf', b, p3), (_Th, gg1, a, s1), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
            - 10 *(
                2 * model.Expand((('Thetakin4W', a, p1, a, p2), (_mf, p1, p2), (_mfa, p2, p1)))
                + model.Expand(((_G, gg1), ('ThetakinWsWf', a, p1), (_Th, gg1, s1, b), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
                + model.Expand(((_G, gg1), ('ThetakinWsWf', a, p2), (_Th, gg1, s1, b), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', a, s1)),MatStruc=['a','a'])
                + model.Expand(((_G, gg1), ('ThetakinWsWf', b, p1), (_Th, gg1, a, s1), (_mf, p1, p2), (_T, gg1, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
                + model.Expand(((_G, gg1), ('ThetakinWsWf', b, p3), (_Th, gg1, a, s1), (_Ta, gg1, p1, p2), (_mf, p2, p3), (_mfa, p3, p1), ('deltatilde', s1, b)),MatStruc=['a','a'])
            )
        )
    res += reskin
    powe[0] = powe[0].subs(BYab, res).subs(BbarYab, 0)
    return powe[0]


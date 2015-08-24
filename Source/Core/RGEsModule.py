#!/usr/bin/env python 
import re as regular  # avoid the double definition of real part of sympy
from RGEsDefinition import *
from ModelsClass import Model, loggingInfo, loggingDebug, loggingCritical
from RGEsmathModule import SimplifyTraces
from Yukawaskin import *
from QuarticCouplingskin import *
from GaugeCouplingskin import *
from FermionMasskin import *
from ScalarMass import *
from Trilinear import *
from ScalarAnomalous import *
from FermionAnomalous import *


def Translate(RGE, model, RunSettings):
    """Do the translation for a given RGE in a general g group into a product of semi simple groups."""
    # One or two loop rge
    LocalRGE = RGEsDefinitions[RGE]['I'] + RGEsDefinitions[RGE]['II'] if RunSettings['Two-Loop'] else \
    RGEsDefinitions[RGE]['I']
    # replace kappa by 1/2 if weyl one otherwise
    LocalRGE = LocalRGE.subs(kappa, Rational(1, 2)) if RunSettings['Weyl'] else LocalRGE.subs(kappa, 1)
    # Separate the different order and store them in Translated
    Translated = {}
    # Check if it is gauge couplings beta functions
    if RGE == 'Gauge-Couplings':
        set_globals(model)
        if model.kinmixing:
            gaugegrouptocalculate = model.NonUGaugeGroups
        else:
            gaugegrouptocalculate = model.GaugeGroups
        for idx, x in enumerate(gaugegrouptocalculate):
            # for idx,x in enumerate(model.GaugeGroups):
            loggingInfo('\t\tGauge-Couplings calculation, Group {} ...'.format(x[0]), verbose=RunSettings['vInfo'])
            # x contains the name and the instance of the Group
            Translated[x[0]] = LocalRGE.expand()
            # separation or the different powers, note the power 0 of g is obtained by putting it to zero
            Translated[x[0]] = [Translated[x[0]].coeff(g ** 2), Translated[x[0]].coeff(g ** 3),
                                Translated[x[0]].coeff(x[1].g ** 4), Translated[x[0]].coeff(g ** 5)]
            # elimination of the powers that do not appear above
            Translated[x[0]] = [[elt * x[1].g ** (ilt + 2), ilt + 2] for ilt, elt in enumerate(Translated[x[0]]) if
                                elt != 0]
            # Loop over the power of g and apply the translation rules
            for ipow, powe in enumerate(Translated[x[0]]):
                if powe[1] == 3:
                    # Eq. 106
                    Translated[x[0]][ipow][0] = ApplyEq106(powe, model, x)
                    loggingInfo('\t\t\t Applying Eq. 106...done', verbose=RunSettings['vInfo'])
                    # Eq. 107
                    Translated[x[0]][ipow][0] = ApplyEq107(powe, model, x)
                    loggingInfo('\t\t\t Applying Eq. 107...done', verbose=RunSettings['vInfo'])
                    # doing the calculation
                    Translated[x[0]][ipow][0] = Translated[x[0]][ipow][0].doit()
                if powe[1] == 5:
                    # Eq. 108
                    Translated[x[0]][ipow][0] = ApplyEq108(powe, model, x)
                    loggingInfo('\t\t\t Applying Eq. 108...done', verbose=RunSettings['vInfo'])
                    # Eq. 109
                    Translated[x[0]][ipow][0] = ApplyEq109(powe, model, x)
                    loggingInfo('\t\t\t Applying Eq. 109...done', verbose=RunSettings['vInfo'])
                    # Eq. 110
                    Translated[x[0]][ipow][0] = ApplyEq110(powe, model, x)
                    loggingInfo('\t\t\t Applying Eq. 110...done', verbose=RunSettings['vInfo'])
                    Translated[x[0]][ipow][0] = Translated[x[0]][ipow][0].doit()
                if powe[1] == 3 and RunSettings['Two-Loop']:
                    # The last step is to calculate the yukawa part
                    Translated[x[0]][ipow][0] = Y4F(powe, model, x)
                    loggingInfo('\t\t\t Calculating the Y4(F) term...done', verbose=RunSettings['vInfo'])
                    Translated[x[0]][ipow][0] = Translated[x[0]][ipow][0].doit()
            Translated[x[0]] = CompileGaugeCouplings(model, Translated, x[0], RunSettings['Weyl'],
                                                     RunSettings['SetGutNorm'])
        # Kinetic mixing terms
        if model.kinmixing:
            loggingInfo('\t\tGauge-Couplings calculation of the abelian sector {} ...'.format(
                ' ,'.join([el[0] for el in model.UGaugeGroups])), verbose=RunSettings['vInfo'])
            # x contains the name and the instance of the Group
            Translated['abelian'] = LocalRGE.expand()
            # separation or the different powers, note the power 0 of g is obtained by putting it to zero
            Translated['abelian'] = [Translated['abelian'].coeff(g ** 2), Translated['abelian'].coeff(g ** 3),
                                     Translated['abelian'].coeff(x[1].g ** 4), Translated['abelian'].coeff(g ** 5)]
            ##elimination of the powers that do not appear above
            Translated['abelian'] = [[elt, ilt + 2] for ilt, elt in enumerate(Translated['abelian']) if elt != 0]
            ##Loop over the power of g and apply the translation rules
            for ipow, powe in enumerate(Translated['abelian']):
                if powe[1] == 3:
                    # Eq. 106
                    Translated['abelian'][ipow][0] = abeliansector1loop(model)
                    loggingInfo('\t\t\t Abelian sector 1-Loop...done', verbose=RunSettings['vInfo'])
                if powe[1] == 5:
                    Translated['abelian'][ipow][0] = abeliansector2loop(model)
                    loggingInfo('\t\t\t Abelian sector 2-Loop...done', verbose=RunSettings['vInfo'])
            Translated['abelian'] = CompileGaugeCouplings(model, Translated, 'abelian', RunSettings['Weyl'],
                                                          RunSettings['SetGutNorm'])
        loggingInfo("\t\t...Gauge-Couplings calculation done.", verbose=RunSettings['vInfo'])
    elif RGE == "Quartic-Couplings":
        set_globalsLbd(model)
        # We loop over the contributions to Calculate and apply the transformation rules
        # Loop over the quartic terms
        if model.LbdToCalculate == {}:
            loggingCritical("WARNING, no Quartic term defined in the lagrangian", verbose=RunSettings['vCritical'])
        else:
            for lbd, listlbd in model.LbdToCalculate.items():
                # we need only one term to calculate the evolution of each lambda let's take the first one in the list
                loggingInfo('\t\tQuartic-Couplings calculation, for parameter {}...\n'.format(lbd),
                            verbose=RunSettings['vInfo'])
                Translated[lbd] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[lbd] = [Translated[lbd].subs(g, 0), Translated[lbd].coeff(g), Translated[lbd].coeff(g ** 2),
                                   Translated[lbd].coeff(g ** 3), Translated[lbd].coeff(g ** 4),
                                   Translated[lbd].coeff(g ** 5), Translated[lbd].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[lbd] = [[elt, ilt] for ilt, elt in enumerate(Translated[lbd]) if elt != 0]
                ######################################
                comb = listlbd[-1]  # select the particles
                ######################################
                for ipow, powe in enumerate(Translated[lbd]):
                    for elemineq in ListEquations['Quartic-Couplings'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[lbd][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['Quartic-Couplings'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[lbd][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[lbd] = CompileQuartic(Translated, lbd, comb, model, RunSettings['Weyl'])
    elif RGE == 'Yukawas':
        set_globalsYuk(model)
        if model.YukToCalculate == {}:
            loggingCritical("WARNING, no Yukawa term defined in the lagrangian", verbose=RunSettings['vCritical'])
        else:
            for y, valy in model.YukToCalculate.items():
                loggingInfo("\t\tParameter {}...".format(y), verbose=RunSettings['vInfo'])
                Translated[y] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[y] = [Translated[y].subs(g, 0), Translated[y].coeff(g), Translated[y].coeff(g ** 2),
                                 Translated[y].coeff(g ** 3), Translated[y].coeff(g ** 4), Translated[y].coeff(g ** 5),
                                 Translated[y].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[y] = [[elt, ilt] for ilt, elt in enumerate(Translated[y]) if elt != 0]
                #####################################################################
                comb = model.YukToCalculate[y][-1]
                sc = [(iel, el) for iel, el in enumerate(comb) if str(el[0]) in model.Scalars]
                assert len(sc) == 1
                comb = [el for iel, el in enumerate(comb) if iel != sc[0][0]]
                comb.insert(0, sc[0][1])
                #####################################################################
                for ipow, powe in enumerate(Translated[y]):
                    for elemineq in ListEquations['Yukawas'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[y][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['Yukawas'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[y][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[y] = CompileYukawas(Translated, y, comb, model, RunSettings['Weyl'],
                                               RunSettings['SetGutNorm'])
                loggingInfo('\t\t Compiling the results for the Yukawa...done', verbose=RunSettings['vInfo'])
    elif RGE == 'ScalarMass':
        set_globalsScM(model)
        if model.ScMToCalculate == {}:
            loggingCritical("WARNING, no Scalar mass term defined in the lagrangian", verbose=RunSettings['vCritical'])
        else:
            for mqt, listmab in model.ScMToCalculate.items():
                #########################################
                # Only one term needed to calculate the evolution lets take the first one The others could be used for implementing some checks
                comb = model.ScMToCalculate[mqt][-1]
                #########################################
                loggingInfo('\t\tScalar mass calculation, for parameter {}...\n'.format(mqt),
                            verbose=RunSettings['vInfo'])
                Translated[mqt] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[mqt] = [Translated[mqt].subs(g, 0), Translated[mqt].coeff(g), Translated[mqt].coeff(g ** 2),
                                   Translated[mqt].coeff(g ** 3), Translated[mqt].coeff(g ** 4),
                                   Translated[mqt].coeff(g ** 5), Translated[mqt].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[mqt] = [[elt, ilt] for ilt, elt in enumerate(Translated[mqt]) if elt != 0]
                for ipow, powe in enumerate(Translated[mqt]):
                    for elemineq in ListEquations['ScalarMass'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[mqt][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['ScalarMass'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[mqt][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[mqt] = CompileScalar(Translated, mqt, comb, model, RunSettings['Weyl'],
                                                RunSettings['SetGutNorm'])
    elif RGE == 'FermionMass':
        set_globalsFM(model)
        if model.FMToCalculate == {}:
            loggingCritical("WARNING, no Fermion mass term defined in the lagrangian", verbose=RunSettings['vCritical'])
        else:
            for mf, valm in model.FMToCalculate.items():
                loggingInfo("\t\tFermion mass calculation started for parameter {}...".format(mf),
                            verbose=RunSettings['vInfo'])
                Translated[mf] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[mf] = [Translated[mf].subs(g, 0), Translated[mf].coeff(g), Translated[mf].coeff(g ** 2),
                                  Translated[mf].coeff(g ** 3), Translated[mf].coeff(g ** 4),
                                  Translated[mf].coeff(g ** 5), Translated[mf].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[mf] = [[elt, ilt] for ilt, elt in enumerate(Translated[mf]) if elt != 0]
                ###############################################################################
                # FOR EACH y i.e. yukawa we have to identify a couple scalar,f1,f2 that gives the yukawa
                comb = model.FMToCalculate[mf][-1]
                ###############################################################################
                for ipow, powe in enumerate(Translated[mf]):
                    for elemineq in ListEquations['FermionMasses'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[mf][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['FermionMasses'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[mf][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[mf] = CompileFM(Translated, mf, comb, model, RunSettings['Weyl'], RunSettings['SetGutNorm'])
                loggingInfo('\t\t Compiling the results for the Yukawa...done', verbose=RunSettings['vInfo'])
    elif RGE == 'Trilinear':
        set_globalsTri(model)
        if model.TriToCalculate == {}:
            loggingCritical("WARNING, no Trilinear term defined in the lagrangian", verbose=RunSettings['vCritical'])
        else:
            for hf, valh in model.TriToCalculate.items():
                loggingInfo('\t\tTrilinear term calculation, for parameter {}...\n'.format(hf),
                            verbose=RunSettings['vInfo'])
                Translated[hf] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[hf] = [Translated[hf].subs(g, 0), Translated[hf].coeff(g), Translated[hf].coeff(g ** 2),
                                  Translated[hf].coeff(g ** 3), Translated[hf].coeff(g ** 4),
                                  Translated[hf].coeff(g ** 5), Translated[hf].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[hf] = [[elt, ilt] for ilt, elt in enumerate(Translated[hf]) if elt != 0]
                # we need only one term to calculate the evolution of each hefg term let's take the first one in the list
                ##############################################################
                comb = model.TriToCalculate[hf][-1]  # Choice only one calculated
                ##############################################################
                for ipow, powe in enumerate(Translated[hf]):
                    for elemineq in ListEquations['Trilinear'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[hf][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['Trilinear'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[hf][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[hf] = CompileTri(Translated, hf, comb, model, RunSettings['Weyl'], RunSettings['SetGutNorm'])
                loggingInfo("\t\t\tCompiling the results for the trilinear terms", verbose=RunSettings['vInfo'])
    elif RGE == 'ScalarAnomalous':
        set_globalsScalarAnomalous(model)
        if model.ScalarAnomalousToCalculate == {}:
            loggingCritical("WARNING, no scalar anomalous dimension to calculate", verbose=RunSettings['vCritical'])
        else:
            for sc, valsc in model.ScalarAnomalousToCalculate.items():
                loggingInfo('\t\tScalar anomalous dimension calculation, for combination {}...\n'.format(sc),
                            verbose=RunSettings['vInfo'])
                Translated[sc] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[sc] = [Translated[sc].subs(g, 0), Translated[sc].coeff(g), Translated[sc].coeff(g ** 2),
                                  Translated[sc].coeff(g ** 3), Translated[sc].coeff(g ** 4),
                                  Translated[sc].coeff(g ** 5), Translated[sc].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[sc] = [[elt, ilt] for ilt, elt in enumerate(Translated[sc]) if elt != 0]
                ##############################################################
                comb = model.ScalarAnomalousToCalculate[sc][-1]
                ##############################################################
                for ipow, powe in enumerate(Translated[sc]):
                    for elemineq in ListEquations['ScalarAnomalous'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[sc][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['ScalarAnomalous'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[sc][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[sc] = CompileScalarAnomalous(Translated, sc, comb, model, RunSettings['Weyl'],
                                                        RunSettings['SetGutNorm'])
                loggingInfo("\t\t\tCompiling the results for the Scalar anomalous dimension terms",
                            verbose=RunSettings['vInfo'])
    elif RGE == 'FermionAnomalous':
        set_globalsFermionAnomalous(model)
        if model.FermionAnomalousToCalculate == {}:
            loggingCritical("WARNING, no fermion anomalous dimension to calculate", verbose=RunSettings['vCritical'])
        else:
            for sc, valsc in model.FermionAnomalousToCalculate.items():
                loggingInfo('\t\tFermion anomalous dimension calculation, for combination {}...\n'.format(sc),
                            verbose=RunSettings['vInfo'])
                Translated[sc] = LocalRGE.expand()
                # separation of the different powers of g
                Translated[sc] = [Translated[sc].subs(g, 0), Translated[sc].coeff(g), Translated[sc].coeff(g ** 2),
                                  Translated[sc].coeff(g ** 3), Translated[sc].coeff(g ** 4),
                                  Translated[sc].coeff(g ** 5), Translated[sc].coeff(g ** 6)]
                # elimination of the component that are equals to zero
                Translated[sc] = [[elt, ilt] for ilt, elt in enumerate(Translated[sc]) if elt != 0]
                ##############################################################
                # comb = model.FermionAnomalousToCalculate[sc][-1]
                comb = [[Symbol('Qbar'), 1, 1], [Symbol('Qbar'), 2, 2]]
                ##############################################################
                for ipow, powe in enumerate(Translated[sc]):
                    for elemineq in ListEquations['FermionAnomalous'][powe[1]]['one-loop']:
                        if elemineq in RunSettings['Skip']:
                            loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                        else:
                            loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                        verbose=RunSettings['vInfo'])
                            Translated[sc][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                            loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                    if RunSettings['Two-Loop']:
                        for elemineq in ListEquations['FermionAnomalous'][powe[1]]['two-loop']:
                            if elemineq in RunSettings['Skip']:
                                loggingCritical('\t\t\t WARNING skipping {} ...'.format(elemineq), verbose=True)
                            else:
                                loggingInfo("\t\t\t Calculating the {} term ...".format(elemineq),
                                            verbose=RunSettings['vInfo'])
                                Translated[sc][ipow][0] = eval("{}(powe,comb,model)".format(elemineq))
                                loggingInfo("\t\t\t\t\t ...done".format(elemineq), verbose=RunSettings['vInfo'])
                Translated[sc] = CompileFermionAnomalous(Translated, sc, comb, model, RunSettings['Weyl'],
                                                         RunSettings['SetGutNorm'])
                loggingInfo("\t\t\tCompiling the results for the Fermion anomalous dimension terms",
                            verbose=RunSettings['vInfo'])
    return Translated

#!/usr/bin/env python
from Contraction import *

try:
    import yaml
    import sys
    from sympy import latex

    sys.path.append('../GroupTheory')
except ImportError:
    exit('error while loading modules in latexRGEs.py')


def writeLatexOutput(RunSettings, FinalRGE, model, ToCalculate):
    try:
        fLatex = open('{}'.format(RunSettings['Results'] + '/' + RunSettings['LatexFile']), 'w')
    except:
        exit('ERROR while opening the latex output file {}'.format(fLatex))
    ModelStr = ''
    ModelStr += 'GaugeGroups :{}\n'.format(model.GaugeGroupsName)
    Particles = copy.deepcopy(model.Particles.values())
    ParticlesNoScalars = [key for key in Particles if not (key.FromCplx)]
    Particles = '\n'.join(['{}: {}'.format(el, el.Qnb) for el in ParticlesNoScalars])
    ModelStr += '\nParticles:\n {}'.format(Particles)
    Potential = '\n\nPotential:'
    for key, val in model.saveSettings['Potential'].items():
        Potential += '\n\n' + key + '\n'
        Potential += '\n'.join(['{}: {}'.format(k, v) for k, v in val.items()])
    ModelStr += Potential
    LatexText = r"""
		\documentclass[12pt,fleqn]{article}
		\pagestyle{empty}
		\usepackage[latin1]{inputenc}
		\usepackage[dvips,landscape,top=1cm,nohead,nofoot]{geometry}
		\usepackage{amsmath}
		\usepackage{amsfonts}
		\usepackage{amssymb}
		\usepackage{scrtime}
		\setlength{\parindent}{0pt}
		\newcommand{\bfrac}[2]{\displaystyle\frac{#1}{#2}}
		\newcommand{\conj}[1]{\ensuremath{#1^*}}
		\newcommand{\transpose}[1]{\ensuremath{#1^t}}
		\newcommand{\trace}[1]{\ensuremath{\mathrm{Tr}#1}}
		\newcommand{\MatM}[1]{\ensuremath{MatM#1}}
		\begin{document}"""
    LatexText += """
		\section{{Summary}}
			\\begin{{flushright}} \\today \\quad \\thistime\end{{flushright}}
			\\begin{{verbatim}}
			model: {0}
			Settings :
			{1}
			\\end{{verbatim}}
			\\newpage
			Model:\\
			\\begin{{verbatim}}
			{2}
			\\end{{verbatim}}
			""".format(RunSettings['Model'], '\n'.join(str(RunSettings).split(','))[1:-1], ModelStr)
    counter = 0
    Translation = {
        'Gauge-Couplings': "writeGaugeCouplings(LatexText,FinalRGE,model,RunSettings)",
        'Quartic-Couplings': "writeQuarticCouplings(LatexText,FinalRGE,model,RunSettings,counter)",
        'Yukawas': "writeYukawas(LatexText,FinalRGE,model,RunSettings,counter)",
        'ScalarMass': "writeScalarMass(LatexText,FinalRGE,model,RunSettings,counter)",
        'FermionMass': "writeFM(LatexText,FinalRGE,model,RunSettings,counter)",
        'Trilinear': "writeTrilinear(LatexText,FinalRGE,model,RunSettings,counter)",
        'ScalarAnomalous': "writeScalarAnomalous(LatexText,FinalRGE,model,RunSettings,counter)",
        'FermionAnomalous': "writeFermionAnomalous(LatexText,FinalRGE,model,RunSettings,counter)"
    }
    for el in ToCalculate:
        #			if model.kinmixing and el =='Gauge-Couplings':
        #				LatexText = writeGaugeCouplingskin(LatexText,FinalRGE,model,RunSettings)
        #			else :
        LatexText = eval(Translation[el])
        counter += 1
    LatexText += '\end{document}'
    # final polishing
    LatexText = LatexText.replace('\\operatorname{transpose}', '\\transpose')
    LatexText = LatexText.replace('\\operatorname{trace}', '\\trace')
    #	LatexText = LatexText.replace('\\operatorname{MatM}','\\MatM')
    fLatex.write(LatexText)
    fLatex.close()


def writeGaugeCouplings(LatexText, FinalRGE, model, RunSettings):
    """writes the GaugeCouplings equation to the latexTexte string"""
    LatexText += '\n\section{Gauge Couplings}'
    LatexText += r"""
		\begin{align*}
			\frac{dg}{dt} = \left.\beta_g\right|_I + \left.\beta_g\right|_{II}\
		\end{align*}"""
    if model.kinmixing:
        toexport = model.NonUGaugeGroups + [('abelian', el, True) for el in flatten(model.UsectorMatrix)]
    else:
        toexport = model.NonUGaugeGroups + [(el[0], el[1], False) for el in model.UGaugeGroups]
    inonu = 0
    for name, g, t in toexport:
        # The Yuk are at two loop but prop to g^3 therefore it ends up in the ExprI contribution -> translate it
        if t:
            ExprI = (FinalRGE[0][name][inonu] * (4 * pi) ** 2).expand().subs(1 / pi, 0)
        else:
            ExprI = (FinalRGE[0][name] * (4 * pi) ** 2).expand().subs(1 / pi, 0)
        # get the dict containing the expression to write on each line
        LtxExprI = SplitEq(ExprI.subs(Symbol('kappa'), Rational(1, 2)).expand())
        gsymbol = latex(g.g) if not (t) else latex(g)
        if RunSettings['Two-Loop']:
            if t:
                ExprII = (FinalRGE[0][name][inonu] * (4 * pi) ** 4).expand().subs(pi ** 2, 0)
            else:
                ExprII = (FinalRGE[0][name] * (4 * pi) ** 4).expand().subs(pi ** 2, 0)
            LtxExprII = SplitEq(ExprII.subs(Symbol('kappa'), Rational(1, 2)).expand(), 200)
            LatexText += r"""
					\subsection{{Evolution of ${0}$}}
					\begin{{align*}}
						\left.(4\pi)^2\beta_{{{0}}}\right|_{{I}} = """.format(gsymbol)
            for line, eq in LtxExprI.items():
                eq = latex(eq)
                Sign = '+' if eq[0] != '-' else ''
                LatexText += r"""
					&{}{}\\""".format(Sign, eq)
            LatexText += r"""\end{align*}
			\begin{align*}"""
            LatexText += r"""
						\left.(4\pi)^4\beta_{{{0}}}\right|_{{II}} = """.format(gsymbol)
            for line, eq in LtxExprII.items():
                eq = latex(eq)
                Sign = '+' if eq[0] != '-' else ''
                LatexText += r"""
				&{}{}\\""".format(Sign, eq)
            LatexText += r"""
					\end{align*}"""
        else:
            LatexText += r"""
					\subsection{{Evolution of ${0}$}}
					\begin{{align*}}
						\left.(4\pi)^2\beta_{{{0}}}\right|_{{I}} = """.format(gsymbol)
            for line, eq in LtxExprI.items():
                eq = latex(eq)
                Sign = '+' if eq[0] != '-' else ''
                LatexText += r"""
					&{}{}\\""".format(Sign, eq)
            LatexText += r"""
					\end{align*}"""
        inonu += 1 if t else 0
    return LatexText


def writeQuarticCouplings(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the Quartic Couplings to the tex string"""
    LatexText += '\n\section{Quartic Coupling}'
    LatexText += r"""
		\begin{align*}
			\frac{d\lambda}{dt} = \left.\beta_{\lambda}\right|_I + \left.\beta_{\lambda}\right|_{II}\
		\end{align*}"""
    index = counter
    for lbda in model.LbdToCalculate.keys():
        Eq = FinalRGE[index][lbda]
        # Get the 1Loop result
        ExprI = ((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2)).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = ((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2)).expand()).subs(pi, 0)
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        if lbda in model.Combination['QuarticTerms']:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, lbda, RunSettings,
                                 combinations=model.Combination['QuarticTerms'][lbda])
        else:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, lbda, RunSettings)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeYukawas(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the Yukawas to the tex string"""
    LatexText += '\n\section{Yukawas}'
    LatexText += r"""
		\begin{align*}
			\frac{dY_a}{dt} = \left.\beta_{Y_a}\right|_I + \left.\beta_{Y_a}\right|_{II}\
		\end{align*}"""
    index = counter
    for y in model.YukToCalculate.keys():
        Eq = FinalRGE[index][y]
        # Reconstruct the external indices
        Indices = ['i{}'.format(ii) for ii in range(len(model.YukToCalculate[y][0][1:]))]
        # Get the 1Loop result
        ExprI = (((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = (((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(pi, 0)
        # Replace the name of the particles with the new indices for readability
        Fermions = [el.args[0] for el in model.YukToCalculate[y][0] if str(el.args[0]) in model.Fermions]
        if len(Fermions) == 2 and Fermions[0] == Fermions[1]:
            ExprI = ExprI.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprI = ExprI.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprI = ExprI.subs(str(Fermions[0]), Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprII = ExprII.subs(str(Fermions[0]), Indices[0])
        else:
            for ipp, pp in enumerate(Fermions):
                ExprI = ExprI.subs(str(pp) + '_f', Indices[ipp])
                ExprI = ExprI.subs(str(pp), Indices[ipp])
                ExprII = ExprII.subs(str(pp) + '_f', Indices[ipp])
                ExprII = ExprII.subs(str(pp), Indices[ipp])
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        if y in model.Combination['Yukawas']:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, y, RunSettings, combinations=model.Combination['Yukawas'][y])
        else:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, y, RunSettings)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeFM(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the FM to the tex string"""
    LatexText += '\n\section{Fermion Mass}'
    LatexText += r"""
		\begin{align*}
			\frac{dm_f}{dt} = \left.\beta_{m_f}\right|_I + \left.\beta_{m_f}\right|_{II}\
		\end{align*}"""
    index = counter
    for m in model.FMToCalculate.keys():
        Eq = FinalRGE[index][m]
        Indices = ['i{}'.format(ii) for ii in range(len(model.FMToCalculate[m][0]))]
        # Get the 1Loop result
        ExprI = (((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = (((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(pi, 0)
        # Split the Equations
        Fermions = [el.args[0] for el in model.FMToCalculate[m][0] if str(el.args[0]) in model.Fermions]
        if len(Fermions) == 2 and Fermions[0] == Fermions[1]:
            ExprI = ExprI.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprI = ExprI.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprII = ExprII.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprI = ExprI.subs(str(Fermions[0]), Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]), Indices[0])
        else:
            for ipp, pp in enumerate(Fermions):
                ExprI = ExprI.subs(str(pp) + '_f', Indices[ipp])
                ExprII = ExprII.subs(str(pp) + '_f', Indices[ipp])
                ExprI = ExprI.subs(str(pp), Indices[ipp])
                ExprII = ExprII.subs(str(pp), Indices[ipp])
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        if m in model.Combination['FermionMasses']:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, m, RunSettings,
                                 combinations=model.Combination['FermionMasses'][m])
        else:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, m, RunSettings)

        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeScalarMass(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the scalar mass to the tex string"""
    LatexText += '\section{Scalar Mass}'
    LatexText += r"""
		\begin{align*}
			\frac{dm_{ab}}{dt} = \left.\beta_{m_{ab}}\right|_I + \left.\beta_{m_{ab}}\right|_{II}\
		\end{align*}"""
    index = counter
    for mab in model.ScMToCalculate.keys():
        Eq = FinalRGE[index][mab]
        # Get the 1Loop result
        ExprI = (((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = (((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(pi, 0)
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        if mab in model.Combination['ScalarMasses']:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, mab, RunSettings,
                                 combinations=model.Combination['ScalarMasses'][mab])
        else:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, mab, RunSettings)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeTrilinear(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the Trilinear terms to the tex string"""
    LatexText += '\section{Trilinear Terms}'
    LatexText += r"""
		\begin{align*}
			\frac{dh_{abc}}{dt} = \left.\beta_{h_{abc}}\right|_I + \left.\beta_{h_{abc}}\right|_{II}\
		\end{align*}"""
    index = counter
    for hfTeX in model.TriToCalculate.keys():
        Eq = FinalRGE[index][hfTeX]
        # Get the 1Loop result
        ExprI = ((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2)).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = ((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2)).expand()).subs(pi, 0)
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        if hfTeX in model.Combination['TrilinearTerms']:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, hfTeX, RunSettings,
                                 combinations=model.Combination['TrilinearTerms'][hfTeX])
        else:
            LatexText = WriteEqs(LatexText, ExprI, ExprII, hfTeX, RunSettings)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeScalarAnomalous(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the scalar anomalous to the tex string"""
    LatexText += '\n\section{Scalar Anomalous dimension}'
    LatexText += r"""
		\begin{align*}
			\frac{d\gamma_{ab}}{dt} = \left.\beta_{\gamma_{ab}}\right|_I + \left.\beta_{\gamma_{ab}}\right|_{II}\
		\end{align*}"""
    index = counter
    for mab in model.ScalarAnomalousToCalculate.keys():
        Eq = FinalRGE[index][mab]
        # Get the 1Loop result
        ExprI = (((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = (((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(pi, 0)
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        #mab = '\\gamma_{{{},{}}}'.format(*[latex(el[0]) for el in model.ScalarAnomalousToCalculate[mab][-1]])
        LatexText = WriteEqs(LatexText, ExprI, ExprII, mab, RunSettings, ScalarAnomalous=True)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def writeFermionAnomalous(LatexText, FinalRGE, model, RunSettings, counter):
    """writes the Fermion anomalous dimension to the tex string"""
    LatexText += '\n\section{Fermion Anomalous dimension}'
    LatexText += r"""
        \begin{align*}
            \frac{d\gamma_{ff^\prime}}{dt} = \left.\beta_{\gamma_{ff^\prime}}\right|_I + \left.\beta_{\gamma_{ff^\prime}}\right|_{II}\\
        \end{align*}
        """
    index = counter
    for y in model.FermionAnomalousToCalculate.keys():
        Eq = FinalRGE[index][y]
        # Reconstruct the external indices
        Indices = ['i{}'.format(ii) for ii in range(2)]
        # Get the 1Loop result
        ExprI = (((Eq * (4 * pi) ** 2).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(1 / pi, 0)
        # Get the 2Loop result
        ExprII = (((Eq * (4 * pi) ** 4).subs(Symbol('kappa'), Rational(1, 2))).expand()).subs(pi, 0)
        # Replace the name of the particles with the new indices for readability
        #Fermions = [el.args[0] for el in model.FermionAnomalousToCalculate[y][0] if str(el.args[0]) in model.Fermions]
        Fermions = [el[0] for el in model.FermionAnomalousToCalculate[y][2] if str(el[0]) in model.Fermions]
        if len(Fermions) == 2 and Fermions[0] == Fermions[1]:
            ExprI = ExprI.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprI = ExprI.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprI = ExprI.subs(str(Fermions[0]), Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]) + '_f', Indices[0])
            ExprII = ExprII.subs(str(Fermions[0]) + '1_f', Indices[1])
            ExprII = ExprII.subs(str(Fermions[0]), Indices[0])
        else:
            for ipp, pp in enumerate(Fermions):
                ExprI = ExprI.subs(str(pp) + '_f', Indices[ipp])
                ExprI = ExprI.subs(str(pp), Indices[ipp])
                ExprII = ExprII.subs(str(pp) + '_f', Indices[ipp])
                ExprII = ExprII.subs(str(pp), Indices[ipp])
        # Split the Equations
        ExprI = SplitEq(ExprI)
        ExprII = SplitEq(ExprII)
        LatexText = WriteEqs(LatexText, ExprI, ExprII, y, RunSettings)
        LatexText = Polishing(LatexText, model)
    else:
        pass
    return LatexText


def WriteEqs(LatexText, ExprI, ExprII, label, RunSettings, combinations=[], ScalarAnomalous=False):
    if RunSettings['Two-Loop']:
        if combinations != []:
            label = latex(Symbol(label)) + '+' + '+'.join(
                [latex(Symbol(el[0]) * el[1]) for el in combinations if el[1] != 0])
            if label[-1] == '+':
                label = latex(label[:-1])
        else:
            if not (ScalarAnomalous):
                label = latex(Symbol(label))
        LatexText += r"""
			\subsection{{Evolution of ${0}$}}
				\begin{{align*}}
					\left.(4\pi)^2\beta_{{{0}}}\right|_{{I}}=""".format(label)
        for line, eq in ExprI.items():
            eq = latex(eq)
            sign = '+' if eq[0] != '-' else ''
            LatexText += r"""
				&{}{}\\""".format(sign, eq)
            if Mod(line, 15) == 0 and line != 0:  # break the block
                LatexText += r"""\end{align*}"""
                LatexText += r"""\begin{align*}
						\phantom{\left.(4\pi)^2\beta_{\lambda}\right|_{I}=}"""
        LatexText += r"""\end{align*}"""
        LatexText += r"""\begin{align*}"""
        LatexText += r"""
			\left.(4\pi)^4\beta_{{{0}}}\right|_{{II}}=""".format(label)
        for line, eq in ExprII.items():
            eq = latex(eq)
            sign = '+' if eq[0] != '-' else ''
            LatexText += r"""
				&{}{}\\""".format(sign, eq)
            if Mod(line, 15) == 0 and line != 0:
                LatexText += r"""\end{align*}"""
                LatexText += r"""\begin{align*}
				\phantom{\left.(4\pi)^2\beta_{\lambda}\right|_{I}=}"""
        LatexText += r"""
			\end{align*}"""
    else:
        if combinations != []:
            label = latex(Symbol(label)) + '+' + '+'.join(
                [latex(Symbol(el[0]) * el[1]) for el in combinations if el[1] != 0])
            if label[-1] == '+':
                label = latex(label[:-1])
        else:
            if not (ScalarAnomalous):
                label = latex(Symbol(label))
        LatexText += r"""
			\subsection{{evolution of ${0}$}}
				\begin{{align*}}
					\left.(4\pi)^2\beta_{{{0}}}\right|_{{I}}=""".format(label)
        for line, eq in ExprI.items():
            eq = latex(eq)
            sign = '+' if eq[0] != '-' else ''
            LatexText += r"""
				&{}{}\\""".format(sign, eq)
            if Mod(line, 15) == 0 and line != 0:
                LatexText += r"""\end{align*}"""
                LatexText += r"""\begin{align*}
				\phantom{\left.(4\pi)^2\beta_{\lambda}\right|_{I}=}"""
        LatexText += r"""
			\end{align*}"""
    return LatexText


def SplitEq(eq, threshold=150):
    """split an equation so that it fits on one line"""
    if type(eq) == Add:
        eq = eq.args
        # count the number of line
        i = 0
        # store the results
        eqout = {}
        eqout[i] = []
        for e in eq:
            tempe = sum(eqout[i]) + e
            # add until we go over the threshold
            if len(latex(tempe)) < threshold:
                eqout[i].append(e)
            else:
                # change the line
                i += 1
                eqout[i] = []
                eqout[i].append(e)
            # recompose the equation
        for k, v in eqout.items():
            eqout[k] = sum(v)
        return eqout
    else:
        eqout = {0: latex(eq)}
        return eqout


def findclosingbracket(lstring, pos):
    ### initialize counter of open brackets
    counter = 0
    while True:
        if lstring[pos] == '}' and counter == 0:
            break
        elif lstring[pos] == '{':
            counter += 1
        elif lstring[pos] == '}':
            counter -= 1
        if pos == len(lstring) - 1 and counter != 0:
            raise Exception('unbalance squarred brackets')
        pos += 1
    return pos


def ExtractSymbs(llist):
    "extract the symbol in between the operatorname brackets"
    pos0 = llist.find('operatorname') + len('operatorname{')
    ll = list(llist)
    pos1 = findclosingbracket(ll, pos0)
    return llist[pos0:pos1]


def Polishing(xpr, model):
    # Simulate the operatorname
    try:
        Yuk = model.ListYukawa + model.ListFM
        OpeName = [ExtractSymbs(latex(model.Classes[el](Symbol(el)))) for el in Yuk]
        # Now we have to substitute '\operatorname{symb} in the whole text
        for iss, symb in enumerate(OpeName):
            xpr = xpr.replace('\operatorname{{{}}}'.format(symb), Yuk[iss])
        xpr = xpr.replace('\overline', '\conj')
    except:
        pass
    return xpr

#!/usr/bin/env python

"""This file contains all the definitions of the RGEs for a simple Lie Group."""

#To get sympy in here
import sys
sys.path.append('../GroupTheory/')
from GroupDefinitions import *
from itertools import permutations
################################
#Symbols used in the definitions
################################
g, F, S, G, kappa, nS = symbols('g F S G kappa,nS')
################
#Gauge couplings 
################
C2, S2, Y4 , Y2= symbols('C2 S2 Y4 Y2',cls = Function)
########################
#Yukawas 1Loop and 2Loop
########################
Yuk, YukGauge, YY2S, YukII, YukGaugeIIg2, YukGaugeIIg4 = symbols('Yuk YukGauge YY2S YukII YukGaugeIIg2 YukGaugeIIg4')
###################
#Quartic invariants
###################
#1Loop
L2abcd,Habcd ,LYabcd = symbols('L2abcd Habcd LYabcd')
#2loop 
L2labcd,Lbar3abcd, Lbar2Yabcd, Hbarlabcd, H2labcd, Hbar2labcd,HYabcd, HbarYabcd,H3abcd = symbols('L2labcd Lbar3abcd Lbar2Yabcd Hbarlabcd H2labcd Hbar2labcd HYabcd HbarYabcd H3abcd')
LSabcd,Lbar2Sabcd,L2gabcd,HSabcd,HFabcd,Y2Flabcd = symbols('LSabcd Lbar2Sabcd L2gabcd HSabcd HFabcd Y2Flabcd ')
#g**4
g4LSabcd, LSSabcd, Alabcd, Abarlabcd, BYabcd, BbarYabcd = symbols(['g4LSabcd', 'LSSabcd', 'Alabcd', 'Abarlabcd', 'BYabcd', 'BbarYabcd'])
#g**6
Atabcd, Agabcd, ASabcd, Aabcd = symbols('Atabcd Agabcd ASabcd Aabcd')
#############
#Fermion Mass
#############
YukFermionMassI,YukFermionMassII,YukGaugeFermionMass,YukGaugeFermionMassIIg2,YukGaugeFermionMassIIg4 = symbols('YukFermionMassI YukFermionMassII YukGaugeFermionMass YukGaugeFermionMassIIg2 YukGaugeFermionMassIIg4')
############
#Scalar Mass
############
meflabef,haefhbef,Hab,LSab,LSabg4,LYab,L2mab,Lbar3ab,Lbar2Yab,Hbarlab,H2mab,Hbar2mab,HYab,HbarYab,H3ab,Lbar2Sab,L2gab,HSab,HFab,Y2Fmab,LSab ,LSSab,Alab,Abarlab,BYab,BbarYab = symbols('meflabef haefhbef Hab LSab LSabg4 LYab L2mab Lbar3ab Lbar2Yab Hbarlab H2mab Hbar2mab HYab HbarYab H3ab Lbar2Sab L2gab HSab HFab Y2Fmab LSab  LSSab Alab Abarlab BYab BbarYab')
#####################
#Trilinear invariants
#####################
#1Loop
L2abc,Habc,LYabc,LSabc = symbols('L2abc Habc LYabc LSabc')
#2 Loop
L2habc,Lbar3abc,Lbar2Yabc,Hbarlmabc,Hbarhabc,H2habc,Hbar2habc,HYabc,HbarYabc,H3abc = symbols('L2habc Lbar3abc Lbar2Yabc Hbarlmabc Hbarhabc H2habc Hbar2habc HYabc HbarYabc H3abc')
Lbar2Sabc,L2gabc,HSabc,HFabc,Y2Fhabc,g4LSabc,LSSabc,Alabc,Abarlabc,BYabc,BbarYabc = symbols('Lbar2Sabc L2gabc HSabc HFabc Y2Fhabc g4LSabc LSSabc Alabc Abarlabc BYabc BbarYabc')
Y2ab,Y2Fab,Hbar2ab,H2ab,ScalarAnoI,ScalarAnoII,xi,L2ab,C2SC2S= symbols('Y2ab Y2Fab Hbar2ab H2ab ScalarAnoI ScalarAnoII xi L2ab C2SC2S')
FermionAnoI,FermionAnoIg2,FermionAnoII,FermionAnoIIg2,FermionAnoIIg4 = symbols('FermionAnoI FermionAnoIg2 FermionAnoII FermionAnoIIg2 FermionAnoIIg4')
####################
#Collect all symbols
####################

ListAllSymbols = {'Yukawas': ['Yuk', 'YukGauge', 'YY2S', 'YukII', 'YukGaugeIIg2', 'YukGaugeIIg4'],
			'QuarticTerms' : ['L2abcd', 'Habcd', 'LYabcd','L2labcd','Lbar3abcd', 'Lbar2Yabcd', 'Hbarlabcd', 'H2labcd' ,'Hbar2labcd', 'HYabcd', 'HbarYabcd', 'H3abcd','LSabcd', 'Lbar2Sabcd', 'L2gabcd', 'HSabcd', 'HFabcd', 'Y2Flabcd','g4LSabcd', 'LSSabcd', 'Alabcd', 'Abarlabcd', 'BYabcd', 'BbarYabcd','Atabcd', 'Agabcd', 'ASabcd', 'Aabcd'],
			'FermionMasses': ['YukFermionMassI','YukFermionMassII', 'YukGaugeFermionMass', 'YukGaugeFermionMassIIg2', 'YukGaugeFermionMassIIg4'],
			'ScalarMasses': ['meflabef', 'haefhbef', 'Hab', 'LSab', 'LSabg4', 'LYab', 'L2mab', 'Lbar3ab', 'Lbar2Yab', 'Hbarlab', 'H2mab', 'Hbar2mab', 'HYab', 'HbarYab', 'H3ab', 'Lbar2Sab', 'L2gab', 'HSab', 'HFab', 'Y2Fmab', 'LSab', 'LSSab', 'Alab', 'Abarlab', 'BYab', 'BbarYab'],
			'TrilinearTerm': ['L2abc', 'Habc', 'LYabc', 'LSabc','L2habc', 'Lbar3abc', 'Lbar2Yabc', 'Hbarlmabc', 'Hbarhabc', 'H2habc', 'Hbar2habc', 'HYabc', 'HbarYabc', 'H3abc','Lbar2Sabc', 'L2gabc', 'HSabc', 'HFabc', 'Y2Fhabc', 'g4LSabc', 'LSSabc', 'Alabc', 'Abarlabc', 'BYabc', 'BbarYabc'],
			'ScalarAnomalous': ['Yab', 'Y2Fab','Hbar2ab', 'H2ab','ScalarAnoI','ScalarAnoII'],
			'FermionAnomalous': ['FermionAnoI','FermionAnoIg2','FermionAnoII','FermionAnoIIg2','FermionAnoIIg4']}


##################
#List of Equations
##################
ListEquations = {}
ListEquations['Quartic-Couplings'] = {0 : {'one-loop': ['CLYabcd','CL2abcd','CHabcd'],
																			'two-loop': ['CL2labcd','CLbar3abcd','CLbar2Yabcd','CHbarlabcd','CHbar2labcd','CH2labcd','CHYabcdHbarYabcdH3abcd']},
																			2 : {'one-loop' : ['CLSabcd'],
																				'two-loop': ['CLbar2Sabcd','CL2gabcd','CY2Flabcd','CHSabcd','CHFabcd']},
																			4 : {'one-loop': ['CAabcd'],
																				'two-loop': ['Cg4LSabcd','CLSSabcd','CAlabcd','CAbarlabcd','CBYabcdBbarYabcd']},
																			6 : {'one-loop': [],
																				'two-loop': ['CtAabcd','CASabcd','CAgabcd']}
																			}
ListEquations['Yukawas'] = {0 : {'one-loop': ['CYuk'],
																'two-loop' : ['CYukII']},
																2 : {'one-loop': ['CYukGauge'],
																	'two-loop' : ['CYukGaugeIIg2']},
																4 : {'one-loop': [],
																	'two-loop': ['CYukGaugeIIg4']}
																}
ListEquations['ScalarMass'] = {0 : {'one-loop': ['Cmeflabef','Chaefhbef','CLYab','CHab'],
																		'two-loop': ['CL2mab','CLbar3ab','CLbar2Yab','CHbarlab','CH2mab','CHbar2mab','CHYabHbarYabH3ab']},
																		2 : {'one-loop': ['CLSab'],
																			'two-loop': ['CLbar2Sab','CL2gab','CY2Fmab','CHSab','CHFab']},
																		4 : {'one-loop': [],
																			'two-loop': ['CLSabg4','CLSSab','CAlab','CAbarlab','CBYabBbarYab']}}
ListEquations['FermionMasses'] = {0 : {'one-loop' :['CYukFermionMassI'],
																			'two-loop': ['CYukFermionMassII']},
																	2 : {'one-loop': ['CYukGaugeFermionMass'],
																		'two-loop': ['CYukGaugeFermionMassIIg2']},
																	4 : {'one-loop': [],
																		'two-loop' : ['CYukGaugeFermionMassIIg4']}}
ListEquations['Trilinear'] = {0 : {'one-loop' :['CL2abc','CHabc','CLYabc'],
																		'two-loop' : ['CL2habc','CLbar3abc','CLbar2Yabc','CHbarlmabc','CHbarhabc','CHbar2habc','CH2habc','CHYabcHbarYabcH3abc']},
																		2 : {'one-loop': ['CLSabc'],
																			'two-loop' : ['CLbar2Sabc','CL2gabc','CY2Fhabc','CHSabc','CHFabc']},
																		4 : {'one-loop': [],
																			'two-loop' : ['Cg4LSabc','CLSSabc','CAlabc','CAbarlabc','CBYabcBbarYabc']
																			}}
ListEquations['ScalarAnomalous'] ={ 0 : {'one-loop': ['CY2ab'],
																					'two-loop': ['CL2ab','CH2ab','CHbar2ab']},
																		2 : {'one-loop': ['CScalarAnoI'],
																					'two-loop': ['CY2Fab']},
																		4 : {'one-loop': [],
																					'two-loop': ['CScalarAnoII','CC2SC2S']}
																		}
ListEquations['FermionAnomalous'] ={ 0 : {'one-loop': ['CFermionAnoI'],
'two-loop': ['CFermionAnoII']},	2 : {'one-loop': ['CFermionAnoIg2'],'two-loop': ['CFermionAnoIIg2']},4 : {'one-loop': [],'two-loop': ['CFermionAnoIIg4']}																		}																	
#################
#RGEs Definitions
#################

RGEsDefinitions = {}
RGEsDefinitions['Gauge-Couplings'] = {
'I':
		-g**3/(4*pi)**2*(
		  Rational(11,3)*C2(G) 
		- Rational(4,3)*kappa*S2(F)
		- Rational(1,6)*S2(S))

		,
'II':
   -g**3/(4*pi)**2*(
   + 2*kappa/(4*pi)**2*Y4(F))

		- g**5/(4*pi)**4*(
		  	Rational(34,3)*C2(G)**2
				- kappa*(4*C2(F)
					+ Rational(20,3)*C2(G))*S2(F)
					-(2*C2(S)
						+ Rational(1,3)*C2(G))*S2(S)
			)
		,
'Ref': 'Eq 30 of 0211440v3'
}
RGEsDefinitions['Quartic-Couplings'] = {
'I':
		1/(4*pi)**2*(
			L2abcd
			- 8*kappa*Habcd
		  + 2*kappa*LYabcd
			- 3*g**2*LSabcd
		  + 3*g**4*Aabcd)
		,
		
'II':
		1/(4*pi)**4*(
			L2labcd
			- Lbar3abcd
			- 4*kappa*Lbar2Yabcd
			+ kappa*(
				8*Hbarlabcd - 3*H2labcd - 2*Hbar2labcd
				)
			+ 4*kappa*(
				HYabcd + 2*HbarYabcd + 2*H3abcd
				)
			+ g**2*(
				2*Lbar2Sabcd - 6*L2gabcd + 4*kappa*(
					HSabcd - HFabcd)
			+ 10*kappa*Y2Flabcd
				)
			- g**4*(
				g4LSabcd - Rational(3,2)*LSSabcd - Rational(5,2)*Alabcd -Rational(1,2)*Abarlabcd + 4*kappa*(BYabcd-10*BbarYabcd)
				)
			+g**6*(Atabcd -Rational(15,2)*ASabcd + 27*Agabcd)
			)
	,
'absnorm': Rational(1,24),
'Ref':'Eq 38 of 0211440v3 and 44'
	}
RGEsDefinitions['Yukawas'] = {
'I':
    1/(4*pi)**2*(
		 Yuk
			-3*g**2*YukGauge)
			,
'II':
	1/(4*pi)**4*(
		YukII
		+ g**2*(YukGaugeIIg2)
		+ g**4*YukGaugeIIg4
		)
	,
'Ref':'Eq 33 and 35 of 0211440v3'
	}
RGEsDefinitions['ScalarMass']={
'I':
  1/(4*pi)**2*(
		 meflabef + haefhbef -4*kappa*Hab -3*g**2*LSab + 2*kappa*LYab)
	,
'II':
  1/(4*pi)**4*(
		 L2mab - Rational(1,2)*Lbar3ab -4*kappa*Lbar2Yab
		 + kappa*(4*Hbarlab - 3*H2mab - 2*Hbar2mab)
     +2*kappa*(HYab +2*HbarYab +2*H3ab)
		 + g**2*(2*Lbar2Sab -6*L2gab + 2* kappa*(HSab-HFab)
			 + 10*kappa*Y2Fmab)
		 -g**4*(
			 LSabg4 - Rational(3,2)*LSSab - Rational(5,2)*Alab - Rational(1,2)*Abarlab + 2*kappa*(BYab-BbarYab))
		 )
	,
'Ref':'Eq 90 of 0211440v3'
}

RGEsDefinitions['FermionMass']={
'I':
    1/(4*pi)**2*(
				YukFermionMassI
				-3*g**2*YukGaugeFermionMass
			)
		,
'II':
     1/(4*pi)**4*(
			 YukFermionMassII
			 + g**2*YukGaugeFermionMassIIg2
			 + g**4*YukGaugeFermionMassIIg4
			 )
		 ,
'Ref':'Eq 62, Eq 63 0f 0211440v3'
}
RGEsDefinitions['Trilinear']={
'I':
   1/(4*pi)**2*(
		 L2abc
		 - 8*kappa*Habc
		 + 2*kappa*LYabc
		 - 3*g**2*LSabc
		 )
	 ,
'II':
   1/(4*pi)**4*(
		 L2habc 
		 - Lbar3abc
		 - 4*kappa*Lbar2Yabc
		 + kappa*(8*Hbarlmabc + 8* Hbarhabc
			 -3*H2habc - 2*Hbar2habc)
		 + 4*kappa*(
			 HYabc + 2*HbarYabc + 2*H3abc)
		 + g**2*(
			 2*Lbar2Sabc -6 * L2gabc + 4* kappa*(HSabc-HFabc)
			 + 10 * kappa* Y2Fhabc)
		 -g**4*(g4LSabc
			 - Rational(3,2)*LSSabc
			 - Rational(5,2)*Alabc
			 - Rational(1,2)*Abarlabc
			 +4*kappa*(BYabc - 10*BbarYabc))
		 )
	 ,
'Ref': 'Eq 64 and 69 of 0211440v3'
}
RGEsDefinitions['ScalarAnomalous'] ={
		'I':
		1/(4*pi)**2*(
			2*kappa*Y2ab 
			-g**2*ScalarAnoI)
		,
		'II':
		1/(4*pi)**4*(
			-g**4*ScalarAnoII
		+ Rational(1,2)*L2ab
		+ Rational(3,2)*g**4*C2SC2S - 3*kappa*H2ab
		- 2 *kappa*Hbar2ab + 10*kappa*g**2*Y2Fab
		)
		,
		'Ref': 'Eq 	22,23 of 0211440v3'
	}	
RGEsDefinitions['FermionAnomalous'] ={
		'I':
		1/(4*pi)**2*(
			FermionAnoI+
		g**2*xi*FermionAnoIg2)
		,
		'II':
		1/(4*pi)**4*(
			FermionAnoII 
			+ g**2*FermionAnoIIg2
			+ g**4*FermionAnoIIg4
		)
		,
		'Ref': 'Eq 	22,23 of 0211440v3'
	}

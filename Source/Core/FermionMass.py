#!/usr/bin/env python 
from RGEsDefinition import *
from ModelsClass import	Model,loggingInfo,loggingDebug,loggingCritical
from RGEsmathModule import DynkinIndex
from RGEsmathModule import DeterminOrdering


###############
#Fermion Masses
###############

def set_globalsFM(model) : 
	"""coppy the globals from RGEsModule into the local gloabls()"""
	for key,val in model.glob.items():
		globals()[key] = val 

def CompileFM(Translated,y,comb,model,Weyl):
	"""Compile a given beta function for the Yukawas i.e. divide by the LHS factor"""
	LH = model.FMToCalculate[y][1]
	FinalBeta = (sum([el[0] for el  in Translated[y]]).doit()/LH).expand()
	if Weyl :
		FinalBeta = FinalBeta.subs(kappa,Rational(1,2)) 
	else :
		FinalBeta = FinalBeta.subs(kappa,1) 
	FinalBeta = FinalBeta.subs(tuple([(el,0) for el in ListAllSymbols['FermionMasses']]))
	FinalBeta  = DeterminOrdering(model,FinalBeta)
	return FinalBeta

######
#1Loop
######

def CYukFermionMassI(powe,comb,model):
	"""Calculate the terms without gauge couplings in Eq 62, first two lines."""
	f1,f2 = comb
	res = (Rational(1,2)*(model.Expand2((([[('Y2Fa',f1,p1),(_mf,p1,f2)],[(_mf,f1,p1),('Y2F',p1,f2)]]))))
			+ 2*model.Expand(((_Y,s1,f1,p1),(_mfa,p1,p2),(_Y,s1,p2,f2)))
				+ kappa*model.Expand2(((_Y,s1,f1,f2),([[(_mfa,p1,p2),(_Y,s1,p2,p1)],[(_mf,p1,p2),(_Ya,s1,p2,p1)]])))
				)
	powe[0] = powe[0].subs(YukFermionMassI,res)
	return powe[0]

#####
#g**2
#####

def CYukGaugeFermionMass(powe,comb,model):
	"""Calculates the part that depends on the gauge couplings line 3 of Eq 62"""
	f1,f2 = comb
	res = model.Expand2(((_G,gg1),([[(_C,gg1,f1),(_mf,f1,f2)],[(_mf,f1,f2),(_C,gg1,f2)]])))
	powe[0] = powe[0].subs(YukGaugeFermionMass,res)
	return powe[0]

######
#2Loop
######

def CYukFermionMassII(powe,comb,model):
	"""Calculates the part that is independent of the gauge couplings in the two loop RGEs for the fermion mass terms Eq.62 8 first lines"""
	f1,f2 = comb	
	res = (
			2*(
				model.Expand(((_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_mf,p2,p3),(_Ya,s1,p3,p4),(_Y,s2,p4,f2)),dotrace=False)#Line 1
			- model.Expand(((_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_mf,p2,p3),(_Ya,s2,p3,p4),(_Y,s1,p4,f2)),dotrace=False)#Line 1
			)
		- model.Expand2(((_Y,s1,f1,p1),([[('Y2F',p1,p2),(_mfa,p2,p3)],[(_mfa,p1,p2),('Y2Fa',p2,p3)]]),(_Y,s1,p3,f2)),dotrace=False)#Line 2

		- Rational(1,8)*(
			model.Expand2((([[(_Y,s1,f1,p1),('Y2F',p1,p2),(_Ya,s1,p2,p3),(_mf,p3,f2)],
				[(_mf,f1,p1),(_Ya,s1,p1,p2),('Y2Fa',p2,p3),(_Y,s1,p3,f2)]])),dotrace=False)#Line 3
			)
		- 2*kappa*(model.Expand2(((_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_Y,s1,p2,f2),([[(_mfa,p3,p4),(_Y,s2,p4,p3)],
			[(_mf,p3,p4),(_Ya,s2,p4,p3)]])))#Line 4 
			)
		- Rational(3,2)*kappa*(model.Expand2((('Yab2S',s1,s2),([[(_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_mf,p2,f2)],
			[(_mf,f1,p1),(_Ya,s2,p1,p2),(_Y,s1,p2,f2)]])),dotrace=False)#Line 5
			)
		- Rational(3,2)*kappa*(model.Expand2(((_Y,s1,f1,f2),([[('Y2F',p1,p2),(_Ya,s1,p2,p3),(_mf,p3,p1)],
			[(_mfa,p1,p2),('Y2Fa',p2,p3),(_Y,s1,p3,p1)]])))#Line 6
			)
		- kappa*(model.Expand2(((_Y,s1,f1,f2),([[(_Y,s2,p1,p2),(_mfa,p2,p3),(_Y,s2,p3,p4),(_Ya,s1,p4,p1)],
			[(_Ya,s2,p1,p2),(_mf,p2,p3),(_Ya,s2,p3,p4),(_Y,s1,p4,p1)]])))#Line 7
			)
		- 2*model.Expand(((_h,s1,s2,s3),(_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_Y,s3,p2,f2)))#Line 8
		+ Rational(1,12)*model.Expand(((_h,s1,s2,s3),(_L,s4,s1,s2,s3),(_Y,s4,f1,f2)))
		)
	powe[0] = powe[0].subs(YukFermionMassII,res)
	return powe[0]


#####
#g**2
#####

def CYukGaugeFermionMassIIg2(powe,comb,model):
	"""Implements Eq 63 lines containing g**2 terms"""
	f1,f2 = comb
	res = (
			3*(model.Expand2(((_G,gg1),([[(_C,gg1,f1),(_Y,s1,f1,p1),(_mfa,p1,p2),(_Y,s1,p2,f2)],
				[(_Y,s1,f1,p1),(_mfa,p1,p2),(_Y,s1,p2,f2),(_C,gg1,f2)]])),dotrace=False)#Line 9
				)
			+ 5*(model.Expand2(((_G,gg1),(_Y,s1,f1,p1),([[(_C,gg1,p1),(_mfa,p1,p2)],[(_mfa,p1,p2),(_C,gg1,p2)]]),(_Y,s1,p2,f2)),dotrace=False)#Line 9
				)
			- Rational(7,4)*(model.Expand2(((_G,gg1),([[(_C,gg1,f1),('Y2Fa',f1,p1),(_mf,p1,f2)],
				[(_mf,f1,p1),('Y2F',p1,f2),(_C,gg1,f2)]])),dotrace=False)#Line 10
				)
			-Rational(1,4)*(model.Expand2(((_G,gg1),([[(_Y,s1,f1,p1),(_C,gg1,p1),(_Ya,s1,p1,p2),(_mf,p2,f2)],
				[(_mf,f1,p1),(_Ya,s1,p1,p2),(_C,gg1,p2),(_Y,s1,p2,f2)]])),dotrace=False)#Line 11
				)
			+ 6*(model.Expand(((_G,gg1),(_Ta,gg1,f1,p1),(_mf,p1,p2),(_Ya,s1,p2,p3),(_Ta,gg1,p3,p4),(_Y,s1,p4,f2)),MatStruc=['a','a'],dotrace=False)#line 12
				+ model.Expand(((_G,gg1),(_Y,s1,f1,p1),(_T,gg1,p1,p2),(_Ya,s1,p2,p3),(_mf,p3,p4),(_T,gg1,p4,f2)),MatStruc=['a','a'],dotrace=False)#Line 12
				)
			+ 5*kappa*(model.Expand2(((_G,gg1),(_Y,s1,f1,f2),(_C,gg1,p1),([[(_mf,p1,p2),(_Ya,s1,p2,p1)],[(_Y,s1,p1,p2),(_mfa,p2,p1)]])))#Line 13
				)
			+ 6*model.Expand(((_G,gg1),(_Cs,gg1,s1),(_Y,s1,f1,p1),(_mfa,p1,p2),(_Y,s1,p2,f2)),dotrace=False)#Line 14

			+	Rational(9,2)*(model.Expand2(((_G,gg1),(_Cs,gg1,s1),([[(_Y,s1,f1,p1),(_Ya,s1,p1,p2),(_mf,p2,f2)],
				[(_mf,f1,p1),(_Ya,s1,p1,p2),(_Y,s1,p2,f2)]])),dotrace=False)#Line 14
					)
			)
	powe[0] = powe[0].subs(YukGaugeFermionMassIIg2,res)
	return powe[0]


#####
#g**4
#####

def CYukGaugeFermionMassIIg4(powe,comb, model):
	"""Calculates the CYukGaugeFermionMAssIIg4"""
	f1,f2 = comb
	res = (-Rational(3,2)*(model.Expand2(((_G,gg1),(_G,gg2),([[(_C,gg1,f1),(_C,gg2,f1),(_mf,f1,f2)],
					[(_mf,f1,f2),(_C,gg1,f2),(_C,gg2,f2)]])),dotrace=False)
				)	
			- Rational(97,6)*(model.Expand2(((_G,gg1),(_G,gg1),(_Cg,gg1),([[(_C,gg1,f1),(_mf,f1,f2)],
					[(_mf,f1,f2),(_C,gg1,f2)]])),dotrace=False)
				)
			+ Rational(10,3)*kappa*(model.Expand2(((_G,gg1),(_G,gg1),(_S,gg1),([[(_C,gg1,f1),(_mf,f1,f2)],
					[(_mf,f1,f2),(_C,gg1,f2)]])),dotrace=False)
				)
			+ Rational(11,12)*(model.Expand2(((_G,gg1),(_G,gg1),(_Ss,gg1),([[(_C,gg1,f1),(_mf,f1,f2)],
					[(_mf,f1,f2),(_C,gg1,f2)]])),dotrace=False)
				)
			)
	powe[0] = powe[0].subs(YukGaugeFermionMassIIg4,res)
	return powe[0]


#!/usr/bin/env python 
import functools,operator,copy
import re as regular #avoid the double definition of real part of sympy 
from RGEsDefinition import *
from ModelsClass import	Model,loggingInfo,loggingDebug,loggingCritical
from RGEsmathModule import DynkinIndex,Epsilon,DeterminOrdering


###########################################################################################################
# Yukawa Transformation rules
###########################################################################################################

def set_globalsYuk(model):
	"""coppy the globals from RGEsModule into the local gloabls()"""
	for key,val in model.glob.items():
		globals()[key] = val 

def CompileYukawas(Translated,y,comb,model,Weyl,GutNorm):
	"""Compile a given beta function for the Yukawas i.e. divide by the LHS factor"""
	LH = model.YukToCalculate[y][1]
	FinalBeta = (sum([el[0] for el  in Translated[y]]).doit()/LH).expand()
	if Weyl :
		FinalBeta = FinalBeta.subs(kappa,Rational(1,2)) 
	else :
		FinalBeta = FinalBeta.subs(kappa,1) 
	if GutNorm:
            FinalBeta = FinalBeta.subs(model.UGaugeGroups[0][1].g,sqrt(Rational(3,5))*model.UGaugeGroups[0][1].g)
	FinalBeta = FinalBeta.subs(tuple([(el,0) for el in ListAllSymbols['Yukawas']]))
	FinalBeta = DeterminOrdering(model,FinalBeta)	
	return FinalBeta

#1loop

def CYuk(powe,comb ,model):
	"""Calculate the term that is independent of the gauge Coupling at one loop"""
	#Eq. 33 
	sc,f1,f2 = comb
	#the first index of comb is the scalar
	res = (Rational(1,2)*(
		  model.Expand2((([[('Y2Fa',f1,p1),(_Y,sc,p1,f2)],
				[(_Y,sc,f1,p1),('Y2F',p1,f2)]])),dotrace=False)
		)
			+ 2*model.Expand(((_Y,s1,f1,p1),(_Ya,sc,p1,p2),(_Y,s1,p2,f2)),dotrace=False)
			+ 2*kappa*(model.Expand(((_Y,s1,f1,f2),('Yab2S',sc,s1)),dotrace=False))
			)
	powe[0] = powe[0].subs(Yuk,res)
	return powe[0]


def CYukGauge(powe,comb,model):
	"""Calculate the term g**2*{C2(F),Ya}= gk(C2_Y * Ya_YZ + Ya_XY *  C2_Y)"""
	sc,f1,f2 = comb
	res = model.Expand2(((_G,gg1),([[(_C,gg1,f1),(_Y,sc,f1,f2)],[(_Y,sc,f1,f2),(_C,gg1,f2)]])),dotrace=False)
	powe[0] = powe[0].subs(YukGauge,res)
	return powe[0]

########
##2Loops
########
def CYukII(powe,comb,model):
	"""Calculates Eq 36 lines 1 to 6"""
	sc,f1,f2 = comb
	res = (2*(model.Expand((('Chain3Y',s2,s1,sc,f1,p1),(_Ya,s2,p1,p2),(_Y,s1,p2,f2)),dotrace=False)
				   - model.Expand((('Chain3Y',s2,s1,sc,f1,p1),(_Ya,s1,p1,p2),(_Y,s2,p2,f2)),dotrace=False)
			  )
				  #Line 1
				 - model.Expand2(((_Y,s1,f1,p1),([[('Y2F',p1,p2),(_Ya,sc,p2,p3)],[(_Ya,sc,p1,p2),('Y2Fa',p2,p3)]]),(_Y,s1,p3,f2)),dotrace=False)#Line 2
				 - Rational(1,8)*(model.Expand2((([[(_Y,s1,f1,p1),('Y2F',p1,p2),(_Ya,s1,p2,p3),(_Y,sc,p3,f2)],
					 	[(_Y,sc,f1,p1),(_Ya,s1,p1,p2),('Y2Fa',p2,p3),(_Y,s1,p3,f2)]])),dotrace=False)#Line 3
						)
				-4*kappa*(model.Expand((('Yab2S',sc,s2),(_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_Y,s1,p2,f2)),dotrace=False))#Line 4
 				-	Rational(3,2)*kappa*(
 						model.Expand2((('Yab2S',s1,s2),([[(_Y,s1,f1,p3),(_Ya,s2,p3,p4),(_Y,sc,p4,f2)],
							[(_Y,sc,f1,p3),(_Ya,s2,p3,p4),(_Y,s1,p4,f2)]])),dotrace=False)#Line 5
 					)
		 	-2*kappa*model.Expand(((_Y,s1,f1,f2),('Hbar2abS',sc,s1)),dotrace=False)#Line 4
				-3*kappa*model.Expand(((_Y,s1,f1,f2),('H2abS',sc,s1)),dotrace=False)#Line 6
				-2*model.Expand(((_L,sc,s1,s2,s3),(_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_Y,s3,p2,f2)),dotrace=False)
				+ Rational(1,2)*model.Expand((('L2abS',sc,s1),(_Y,s1,f1,f2)),dotrace=False)
			)
	powe[0] = powe[0].subs(YukII,res)
	return powe[0]

def CYukGaugeIIg2(powe,comb,model):
	"""Calculates Eq 36 Lines 7-12"""
	sc,f1,f2 = comb
	res = (
			3*(model.Expand2(((_G,gg1),([[(_C,gg1,f1),(_Y,s1,f1,p1),(_Ya,sc,p1,p2),(_Y,s1,p2,f2)],
				[(_Y,s1,f1,p1),(_Ya,sc,p1,p2),(_Y,s1,p2,f2),(_C,gg1,f2)]])),dotrace=False)
			)#Line 7
		+ 5*(model.Expand2(((_G,gg1),(_Y,s1,f1,p1),([[(_C,gg1,p1),(_Ya,sc,p1,p2)],[(_Ya,sc,p1,p2),(_C,gg1,p2)]]),(_Y,s1,p2,f2)),dotrace=False)
			)#Line 7
		- Rational(7,4)*(model.Expand2(((_G,gg1),([[(_C,gg1,f1),('Y2Fa',f1,p1),(_Y,sc,p1,f2)],[(_Y,sc,f1,p1),('Y2F',p1,f2),(_C,gg1,f2)]])),dotrace=False)
			)#Line 8
		-Rational(1,4)*(model.Expand2(((_G,gg1),([[(_Y,s1,f1,p1),(_C,gg1,p1),(_Ya,s1,p1,p2),(_Y,sc,p2,f2)],
				[(_Y,sc,f1,p1),(_Ya,s1,p1,p2),(_C,gg1,p2),(_Y,s1,p2,f2)]])),dotrace=False)
			)#Line 9
		+ 10*kappa*model.Expand(((_G,gg1),(_Y,s1,f1,f2),('Y2FabS',gg1,sc,s1)),dotrace=False)#Line 10
		+ 6*(model.Expand(((_G,gg1),(_Cs,gg1,s1),(_Y,s1,f1,p1),(_Ya,sc,p1,p2),(_Y,s1,p2,f2)),dotrace=False)#Line 11
		   -2*model.Expand(((_G,gg1),(_Cs,gg1,sc),(_Y,s1,f1,p1),(_Ya,sc,p1,p2),(_Y,s1,p2,f2)),dotrace=False)#Line 11
			)
		+ Rational(9,2)*(model.Expand2(((_G,gg1),(_Cs,gg1,s1),([[(_Y,s1,f1,p1),(_Ya,s1,p1,p2),(_Y,sc,p2,f2)],
				[(_Y,sc,f1,p1),(_Ya,s1,p1,p2),(_Y,s1,p2,f2)]])),dotrace=False)
			)
		+ 6*(model.Expand(((_G,gg1),(_Ta,gg1,f1,p1),(_Y,sc,p1,p2),(_Ya,s1,p2,p3),(_Ta,gg1,p3,p4),(_Y,s1,p4,f2)),MatStruc=['A','A'])
		 			+ model.Expand(((_G,gg1),(_Y,s1,f1,p1),(_T,gg1,p1,p2),(_Ya,s1,p2,p3),(_Y,sc,p3,p4),(_T,gg1,p4,f2)),MatStruc=['A','A'])
			)
	)
	powe[0] = powe[0].subs(YukGaugeIIg2,res)
	return powe[0]

def CYukGaugeIIg4(powe,comb,model):
	"""Calculates Eq 36 Lines 13-End"""
	sc,f1,f2 = comb
	res = (
			-Rational(3,2)*(model.Expand2(((_G,gg1),(_G,gg2),([[(_C,gg1,f1),(_C,gg2,f1),(_Y,sc,f1,f2)],
				[(_Y,sc,f1,f2),(_C,gg1,f2),(_C,gg2,f2)]])),dotrace=False)#Line 13
				)
			+ 6*(model.Expand2(((_G,gg1),(_G,gg2),(_Cs,gg1,sc),([[(_C,gg2,f1),(_Y,sc,f1,f2)],[(_Y,sc,f1,f2),(_C,gg2,f2)]])),dotrace=False)#Line 13
				)
			- Rational(97,6)*(model.Expand2(((_G,gg1),(_G,gg1),(_Cg,gg1),([[(_C,gg1,f1),(_Y,sc,f1,f2)],
					[(_Y,sc,f1,f2),(_C,gg1,f2)]])),dotrace=False)#Line 14
				)
			+ Rational(10,3)*kappa*(model.Expand2(((_G,gg1),(_G,gg1),(_S,gg1),([[(_C,gg1,f1),(_Y,sc,f1,f2)],
					[(_Y,sc,f1,f2),(_C,gg1,f2)]])),dotrace=False)#Line 14
				)
			+ Rational(11,12)*(model.Expand2(((_G,gg1),(_G,gg1),(_Ss,gg1),([[(_C,gg1,f1),(_Y,sc,f1,f2)],[(_Y,sc,f1,f2),(_C,gg1,f2)]])),dotrace=False)#Line 15
				)
			- Rational(21,2)*(model.Expand(((_G,gg1),(_G,gg2),(_Cs,gg1,sc),(_Cs,gg2,sc),(_Y,sc,f1,f2)),dotrace=False))#Line 15
			+ Rational(49,4)*model.Expand(((_G,gg1),(_G,gg1),(_Cs,gg1,sc),(_Cg,gg1),(_Y,sc,f1,f2)),dotrace=False)#Line 16
			- 2*kappa*model.Expand(((_G,gg1),(_G,gg1),(_Cs,gg1,sc),(_S,gg1),(_Y,sc,f1,f2)),dotrace=False)#Line 16
			- Rational(1,4)*model.Expand(((_G,gg1),(_G,gg1),(_Cs,gg1,sc),(_Ss,gg1),(_Y,sc,f1,f2)),dotrace=False)#Line 16
			)
	powe[0] = powe[0].subs(YukGaugeIIg4,res)
	return powe[0]

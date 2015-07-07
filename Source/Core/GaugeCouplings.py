#!/usr/bin/env python 
import functools,operator,copy
import re as regular #avoid the double definition of real part of sympy 
from RGEsDefinition import *
from ModelsClass import	Model,loggingInfo,loggingDebug,loggingCritical
from RGEsmathModule import DynkinIndex,DynkinCasimir,DeterminOrdering
################################################################################
#Gauge Couplings ransformation rules 
################################################################################
#update global variables 

#powe 3
def CompileGaugeCouplings(model,Translated,x,Weyl): 
	FinalBeta = sum([el[0] for el in Translated[x]]).expand()
	if Weyl :
		FinalBeta = FinalBeta.subs(kappa,Rational(1,2))
	else :
		FinalBeta = FinalBeta.subs(kappa,1)
	FinalBeta = DeterminOrdering(model,FinalBeta)
	return FinalBeta

def set_globals(model):
	"""coppy the globals from RGEsModule into the local gloabls()"""
	for key,val in model.glob.items():
		globals()[key] = val 

def ApplyEq106(powe, Mod,x) :
	"""Apply the transformation rule defined by Eq. 106"""
	assert powe[1] == 3
	if not(x[1].U) : 
		return powe[0].subs([(C2(G),x[1].C2('G'))]) if powe[1] == 3 else powe[0]
	else :
		"There is no Casimir for the adjoint for the U1 group"
		return powe[0].subs([(C2(G),0)]) if powe[1] == 3 else powe[0]

def ApplyEq107(powe, Mod, x) :
	"""Apply the transformation rule defined by Eq. 107"""
	assert powe[1] == 3
	out = powe[0].subs([(S2(F),DynkinIndex(Mod, x,'Fermions')),(S2(S),DynkinIndex(Mod, x,'Scalars'))])
	return out

#power 5
def ApplyEq108(powe, Mod, x) :
	"""Apply the transformation rule defined by Eq. 108."""
	assert powe[1] == 5
	if not(x[1].U) :
		loggingDebug('Applying equation 108 on : {}'.format(powe))
		return powe[0].subs(C2(G)**2,(x[1].C2('G'))**2) 
	else :
		return powe[0].subs(C2(G),0)

def ApplyEq109(powe, Mod, x) :
	"""Apply the transformation rule defined by Eq. 109."""
	assert powe[1] == 5
	if not(x[1].U) :
		return powe[0].subs([(C2(G)*S2(F),x[1].C2('G')*DynkinIndex(Mod,x,'Fermions')),(C2(G)*S2(S),x[1].C2('G')*DynkinIndex(Mod,x,'Scalars'))])
		loggingDebug('Applying equation 109 on : {}'.format(powe))
	else :
		return powe[0].subs([(C2(G)*S2(F),0),(C2(G)*S2(S),0)])

def ApplyEq110(powe, Mod, x) :
	"""Apply the transformation rule defined by Eq. 110."""
	assert powe[1] == 5
	return powe[0].expand().subs([(x[1].g**5*C2(F)*S2(F),DynkinCasimir(Mod, x, 'Fermions')),(x[1].g**5*C2(S)*S2(S),DynkinCasimir(Mod, x, 'Scalars'))])


def Y4F(powe, model, x):
	"""Calculates the yukawa term in the gauge coupling, i.e. Y4(F)."""
	assert powe[1] == 3
	# Expand the product 1/d*C2(i)*Y(z,i,j)*Ydagger(z,j,i) z among higgs fields and i,j among fermions
	res = Rational(1,x[1].d)*model.Expand(((_C,x[0],p1),(_Y,s1,p1,p2),(_Ya,s1,p2,p1)))
	return powe[0].expand().subs(Y4(F),res)



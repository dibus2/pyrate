#!/usr/bin/env python 
import functools,operator,copy
import re as regular #avoid the double definition of real part of sympy 
from RGEsDefinition import *
from ModelsClass import	Model,loggingInfo,loggingDebug,loggingCritical
from RGEsmathModule import DynkinIndex,DynkinCasimir,multiplicity,SimplifyTraces,DeterminOrdering,Tr
import time
import pudb


##################################################################################
#Quartic coupling functions 
##################################################################################

def set_globalsLbd(model):
	"""coppy the globals from RGEsModule into the local gloabls()"""
	for key,val in model.glob.items():
		globals()[key] = val 

def CompileQuartic(Translated,lbd,comb,model,Weyl):
	LH = model.LbdToCalculate[lbd][1]
	FinalBeta = (sum([el[0] for el in Translated[lbd]]).doit()/LH).expand()
	if Weyl : 
		FinalBeta = FinalBeta.subs(kappa,Rational(1,2)) 
	else :
		FinalBeta = FinalBeta.subs(kappa,1) 
	#Set all the terms that have been skipped to zero
	FinalBeta = FinalBeta.subs(tuple([(el,0) for el in ListAllSymbols['QuarticTerms']]))
	FinalBeta = DeterminOrdering(model,FinalBeta)	
	return FinalBeta
#power 0 

def CLYabcd(powe,comb,model):
	"""Calculate the LYabcd term."""
	assert powe[1] == 0
	sc1,sc2,sc3,sc4 = comb
	res = [model.Expand((('Yab2S',sc,sc),(_L,sc1,sc2,sc3,sc4)))
			for sc in comb]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(LYabcd,res)
	return powe[0]

def CL2abcd(powe,comb,model):
	"""Calculate the L2abcd term"""
	assert powe[1] == 0
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,8)*model.Expand(((_L,a,b,s1,s2),(_L,s1,s2,c,d)))
				for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
				]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(L2abcd,res)
	return powe[0]


def CHabcd(powe,comb,model):
	"""Calculates the Habcd term Eq 40"""
	assert powe[1] == 0
	sc1,sc2,sc3,sc4 = comb
	res = model.Expand((('Habcd',sc1,sc2,sc3,sc4)))
	powe[0] = powe[0].subs(Habcd,res)
	return powe[0]

#power 2 
def CLSabcd(powe,comb,model):
	"""Calculate the term LSabcd Eq 42.
	This gets modified by the kin mixing according to Eq. (17)"""
	sc1,sc2,sc3,sc4 = comb
	res = [model.Expand(((_G,gg1),(_Cs,gg1,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4]
			]
	res = sum(res)
	reskin = 0 
	if model.kinmixing :
			reskin = [model.Expand2(((_Ckins,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4]
			]
			reskin = sum(reskin)
	res = res + reskin
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(LSabcd,res)
	return powe[0]

#power 4 
def CAabcd(powe,comb, model):
	"""Calculate the Aabcd term Eq 43 + replacement"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,4)*(
		model.Expand(((_G,gg1),(_G,gg2),('Theta4',gg1,gg2,a,b,c,d)))
			)
		for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
		]
	res = sum(res)
	reskin = 0 
	if model.kinmixing : 
			reskin = [(
					#Terms Lambda*C
					model.Expand((('Theta2g2kin',a,c,s1,s2),('Thetakin',s1,b,s2,d)))
					+ model.Expand((('Theta2g2kin',s1,s2,b,d),('Thetakin',a,s1,c,s2)))
					+ model.Expand((('Theta2g2kin',a,s2,s1,d),('Thetakin',s1,b,c,s2)))
					+ model.Expand((('Theta2g2kin',s1,c,b,s2),('Thetakin',a,s1,s2,d)))
					#Terms C*C
					+ model.Expand((('Thetakin',a,s1,c,s2),('Thetakin',s1,b,s2,d)))
					+ model.Expand((('Thetakin',a,s1,s2,d),('Thetakin',s1,b,c,s2)))
					)
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
			reskin = sum(reskin)
	res = res + reskin 	
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Aabcd,res)
	return powe[0]


#####
#2Loop
#####

##########
#Invariant
##########

#########
#labcd**3
#########


def CL2labcd(powe,comb, model):
	"""Calculates the term sum_L2(i)*labcd Eq 24 and 44"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,2)*model.Expand((('L2abS',a,a),(_L,sc1,sc2,sc3,sc4)),dotrace=False)
			for a in [sc1,sc2,sc3,sc4]
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(L2labcd,res)
	return powe[0]

def CLbar3abcd(powe,comb, model):
	"""Calculates the term Lbar3abcd Eq 44 lign 1"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,4)*model.Expand(((_L,a,b,s1,s2),(_L,c,s1,s3,s4),(_L,d,s2,s3,s4)))
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Lbar3abcd,res)
	return powe[0]

#########
#labcd**2
#########

def CLbar2Yabcd(powe,comb,model):
	"""Calculates the term CLbar2Yabcd Eq 46."""
	sc1,sc2,sc3,sc4 = comb
 	res =	[Rational(1,8)*(model.Expand((('Yab2S',s1,s2),(_L,a,b,s3,s1),(_L,c,d,s3,s2))))
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Lbar2Yabcd,res)
	return powe[0]

######
#labcd
######

def CHbarlabcd(powe,comb, model):
	"""Calculates the term Hbarlabcd Eq 47"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,8)*(
			#model.Expand(((_L,a,b,s1,s2),('Chain4Y',c,s1,d,s2,p1,p1)))
			model.Expand(((_L,a,b,s1,s2),('Chain3Y',c,s1,d,p1,p2),(_Ya,s2,p2,p1)))
			#+model.Expand(((_L,a,b,s1,s2),('Chain4Ya',c,s1,d,s2,p1,p1)))
			+model.Expand(((_L,a,b,s1,s2),(_Ya,c,p1,p2),('Chain3Y',s1,d,s2,p2,p1)))
		)
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Hbarlabcd,res)
	return powe[0]

def CHbar2labcd(powe,comb,model):
	"""Calculates the term Hbar2labcd Eq 25"""
	sc1,sc2,sc3,sc4 = comb
	res =[ model.Expand((('Hbar2abS',i,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4]
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Hbar2labcd,res)
	return powe[0]

def CH2labcd(powe,comb, model):
	"""Calculates the term H2labcd Eq 25"""
	sc1,sc2,sc3,sc4 = comb
	res =[ model.Expand((('H2abS',i,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4]
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(H2labcd,res)
	return powe[0]



####
#Yuk
####

def CHYabcdHbarYabcdH3abcd(powe,comb, model):
	"""Calculates the terms HYabcd HbarYabcd and H3abcd to speed up the calculation Eqs 48-49-50"""
	sc1,sc2,sc3,sc4 = comb
	calculated = {}
	tot =[]
	for perm in list(permutations([sc1,sc2,sc3,sc4],4)):
		key = tuple(flatten(perm))
		if key in calculated : 
			res = calculated[key]
			tot.append(res)
		else :
			a,b,c,d = perm 
			resHY = (model.Expand((('Y2F',p1,p2),(_Ya,a,p2,p3),('Chain3Y',b,c,d,p3,p1))))#HY 
			resHYbar =  2*Rational(1,2)*(#HbarY
							model.Expand(((_Y,s1,p1,p2),(_Ya,a,p2,p3),('Chain3Y',s1,b,c,p3,p4),(_Ya,d,p4,p1)))
							)
			#resHYbarad =  2*Rational(1,2)*(#HbarY
			#				model.Expand(((_Y,d,p1,p2),(_Ya,c,p2,p3),('Chain3Y',b,s1,a,p3,p4),(_Ya,s1,p4,p1)))
			#				)
			resH3 = 2*Rational(1,2)*(model.Expand((('Chain3Y',a,b,s1,p1,p2),(_Ya,c,p2,p3),(_Y,d,p3,p4),(_Ya,s1,p4,p1)))
						)
			resH3 = resH3.doit() if resH3 != 0 else 0
			resHY = resHY.doit() if resHY != 0 else 0
			resHYbar = resHYbar.doit() if resHYbar != 0 else 0 
			resHYbarad = resHYbar.adjoint() if resHYbar != 0 else 0 
			resHYbarad = resHYbarad.doit() if resHYbarad != 0 else 0
			res = SimplifyTraces(resHYbar + resHYbarad + resHY + resH3,model)
			if res != 0 : 
				tot.append(res)
			calculated[key] = res
	res = sum(tot)
	powe[0] = powe[0].subs(HYabcd,res).subs(HbarYabcd,0).subs(H3abcd,0)
	return powe[0]

#####
#g**2
#####

def CLbar2Sabcd(powe,comb,model):
	"""Calculates the term CLbar2Yabcd Eq 51."""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,8)*model.Expand(((_G,gg1),(_Cs,gg1,s1),(_L,a,b,s2,s1),(_L,c,d,s2,s1)))
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res) .doit()
	powe[0] = powe[0].subs(Lbar2Sabcd,res)
	return powe[0]

def CL2gabcd(powe,comb,model):
	"""Calculates the term L2gabcd Eq52 """
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,8)*model.Expand(((_G,gg1),(_L,a,b,s1,s2),(_L,c,d,s3,s4),(_Th,gg1,s1,s3),(_Th,gg1,s2,s4)),MatStruc=['A','A'])
			for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(L2gabcd,res)
	return powe[0]

def CY2Flabcd(powe,comb,model):
	"""Calculates the term Y2Flabcd Eq27 """
	sc1,sc2,sc3,sc4 = comb
	res = [model.Expand(((_G,gg1),('Y2FabS',gg1,i,i),(_L,sc1,sc2,sc3,sc4)))
		for i in [sc1,sc2,sc3,sc4]
		]
	res = sum(res)
	if res != 0 :
		res = res.doit()

	powe[0] = powe[0].subs(Y2Flabcd,res)
	return powe[0]

def CHSabcd(powe,comb, model):
	"""Calculates the HSabcd term Eq 53"""
	sc1,sc2,sc3,sc4 = comb
	res = [model.Expand(((_G,gg1),(_Cs,gg1,i),('Habcd',sc1,sc2,sc3,sc4)))
		for i in [sc1,sc2,sc3,sc4]
		]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(HSabcd,res)
	return powe[0]

def CHFabcd(powe,comb, model):
	"""Calculates the term HFabcd Eq 54"""
	sc1,sc2,sc3,sc4 = comb
	res =[model.Expand2(((_G,gg1),([[(_C,gg1,p1),(_Y,a,p1,p2)],[(_Y,a,p1,p2),(_C,gg1,p2)]]),('Chain3Ya',b,c,d,p2,p1)))
			for  (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(HFabcd,res)
	return powe[0]

#####
#g**4
#####


def Cg4LSabcd(powe,comb,model):
	"""Calculate the term LSabcd Eq 44 line 6"""
	sc1,sc2,sc3,sc4 = comb
	res =[ Rational(35,3)*model.Expand(((_G,gg1),(_G,gg1),(_Cg,gg1),(_Cs,gg1,i),(_L,sc1,sc2,sc3,sc4)))
			 -Rational(10,3)*kappa*model.Expand(((_G,gg1),(_G,gg1),(_S,gg1),(_Cs,gg1,i),(_L,sc1,sc2,sc3,sc4)))
			 -Rational(11,12)*model.Expand(((_G,gg1),(_G,gg1),(_Ss,gg1),(_Cs,gg1,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4 ]
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(g4LSabcd,res)
	return powe[0]


def CLSSabcd(powe,comb,model):
	"""Calculate the term LSSabcd Eq 55"""
	sc1,sc2,sc3,sc4 = comb
	res = [model.Expand(((_G,gg1),(_G,gg2),(_Cs,gg1,i),(_Cs,gg2,i),(_L,sc1,sc2,sc3,sc4)))
			for i in [sc1,sc2,sc3,sc4 ]
			]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(LSSabcd,res)
	return powe[0]

def CAlabcd(powe,comb,model):
	"""Calculates the tern Alabcd Eq 56"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,2)*(
		model.Expand(((_G,gg1),(_G,gg2),(_L,a,b,s1,s2),('Theta4',gg1,gg2,s1,s2,c,d)))
	)
		for  (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
		]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Alabcd,res)
	return powe[0]
		
def CAbarlabcd(powe,comb,model):
	"""Calculates the tern Albarabcd Eq 56"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,2)*(
		model.Expand(((_G,gg1),(_G,gg2),(_L,a,b,s1,s2),('Theta4',gg1,gg2,c,s1,d,s2)))
	)
		for  (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
		]
	res = sum(res)
	if res != 0 :
		res = res.doit()
	powe[0] = powe[0].subs(Abarlabcd,res)
	return powe[0]

									
def CBYabcdBbarYabcd(powe,comb,model):
	"""Calculates the term BYabcd Eq 58"""
	sc1,sc2,sc3,sc4 = comb
	tot = []
	calculated = {}
	ts = time.time()
	By,Bybar = symbols('By Bybar')
	for perm in list(permutations([sc1,sc2,sc3,sc4],4)):
		key = tuple(flatten(perm))
		if not(key in calculated): 
			a,b,c,d = perm
			res = Rational(1,4)*(
			model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_Ta,gg4,p2,p3),(_Y,c,p3,p4),(_Ya,d,p4,p1)),MatStruc=[['a','b','a','b'],['b','a','a','b']])
			+ model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Y,c,p1,p2),(_T,gg3,p2,p3),(_T,gg4,p3,p4),(_Ya,d,p4,p1)),MatStruc=[['a','b','a','b'],['b','a','a','b']])
			-10*model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_Y,c,p2,p3),(_T,gg4,p3,p4),(_Ya,d,p4,p1)),MatStruc=[['a','b','a','b'],['b','a','a','b']])
			)
			if res != 0 : 
				tot.append(res.doit())
			calculated[key] = res.doit()
		else : 
			tot.append(calculated[key])
	res = sum(tot)
	res = res.doit()
	tf = time.time()
	powe[0] = powe[0].subs(BYabcd,res).subs(BbarYabcd,0)
	return powe[0]

#def CBYabcdBbarYabcd(powe,comb,model):
#	"""Calculates the term BYabcd Eq 58"""
#	sc1,sc2,sc3,sc4 = comb
#	tot = []
#	calculated = {}
#	By,Bybar = symbols('By Bybar')
#	for perm in list(permutations([sc1,sc2,sc3,sc4],4)):
#		key = tuple(flatten(perm))
#		if not(key in calculated): 
#			a,b,c,d = perm
#			res = Rational(1,4)*(
#			model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg1,p1,p2),(_Ta,gg2,p2,p3),(_Y,c,p3,p4),(_Ya,d,p4,p1)),MatStruc=['a','b','a','b'])
#			+ model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg2,p1,p2),(_Ta,gg1,p2,p3),(_Y,c,p3,p4),(_Ya,d,p4,p1)),MatStruc=['b','a','a','b'])
#			+ model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Y,c,p1,p2),(_T,gg1,p2,p3),(_T,gg2,p3,p4),(_Ya,d,p4,p1)),MatStruc=['a','b','a','b'])
#			+ model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Y,c,p1,p2),(_T,gg2,p2,p3),(_T,gg1,p3,p4),(_Ya,d,p4,p1)),MatStruc=['b','a','a','b'])
#			-10*model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg1,p1,p2),(_Y,c,p2,p3),(_T,gg2,p3,p4),(_Ya,d,p4,p1)),MatStruc=['a','b','a','b'])
#			-10*model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg2,p1,p2),(_Y,c,p2,p3),(_T,gg1,p3,p4),(_Ya,d,p4,p1)),MatStruc=['b','a','a','b'])
#			)
#			if res != 0 : 
#				tot.append(res.doit())
#			calculated[key] = res.doit()
#		else : 
#			tot.append(calculated[key])
#	res = sum(tot)
#	powe[0] = powe[0].subs(BYabcd,res).subs(BbarYabcd,0)
#	return powe[0]


#######
#g**6
#######

def CtAabcd(powe,comb,model):
	"""Calculates the term tAabcd Eq (44) last line"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(161,6)*Rational(1,4)*(
			model.Expand(((_G,gg1),(_G,gg1),(_G,gg2),(_Cg,gg1),('Theta4',gg1,gg2,a,b,c,d)))
		)
		- Rational(32,3)*kappa*Rational(1,4)*(
			model.Expand(((_G,gg1),(_G,gg1),(_G,gg2),(_S,gg1),('Theta4',gg1,gg2,a,b,c,d)))
		)
		-Rational(7,3)*Rational(1,4)*(
			model.Expand(((_G,gg1),(_G,gg1),(_G,gg2),(_Ss,gg1),('Theta4',gg1,gg2,a,b,c,d)))
		)
		for  (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
		]
	res = sum(res)
	powe[0] = powe[0].subs(Atabcd,res)
	return powe[0]

def CASabcd(powe,comb,model):
	"""Calculate the term ASabcd Eq 60"""
	sc1,sc2,sc3,sc4 = comb
	res = [Rational(1,4)*(
			model.Expand(((_G,gg1),(_Cs,gg1,ii),(_G,gg2),(_G,gg3),('Theta4',gg2,gg3,a,b,c,d)))
		)
		for ii in [sc1,sc2,sc3,sc4]
		for  (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
		]
	res = sum(res)
	powe[0] = powe[0].subs(ASabcd,res)
	return powe[0]


def CAgabcd(powe,comb,model):
	"""Calculate the term Agabcd Eq 61"""
	sc1,sc2,sc3,sc4 = comb
	tot = []
	calculated = {}
	for perm in list(permutations([sc1,sc2,sc3,sc4],4)):
		key = tuple(flatten(perm))
		if not(key in calculated):
			a,b,c,d = perm
			res = Rational(1,8)*(
				model.Expand(((_G,gg1),(_G,gg1),(_G,gg1),('ffabc',gg1),(_Th,gg1,a,s1),(_Th,gg1,s1,b),(_Th,gg1,c,s2),(_Th,gg1,s2,d)),MatStruc=['a','c','e','b','d','e','a','b','c','d'])
			+ model.Expand(((_G,gg1),(_G,gg1),(_G,gg1),('ffabc',gg1),(_Th,gg1,a,s1),(_Th,gg1,s1,b),(_Th,gg1,c,s2),(_Th,gg1,s2,d)),MatStruc=['a','c','e','b','d','e','a','b','d','c'])
			+ model.Expand(((_G,gg1),(_G,gg1),(_G,gg1),('ffabc',gg1),(_Th,gg1,a,s1),(_Th,gg1,s1,b),(_Th,gg1,c,s2),(_Th,gg1,s2,d)),MatStruc=['a','c','e','b','d','e','b','a','d','c'])
			+ model.Expand(((_G,gg1),(_G,gg1),(_G,gg1),('ffabc',gg1),(_Th,gg1,a,s1),(_Th,gg1,s1,b),(_Th,gg1,c,s2),(_Th,gg1,s2,d)),MatStruc=['a','c','e','b','d','e','b','a','c','d'])
					)
			if res != 0 : 
				tot.append(res)
			calculated[key] = res
		else  :
			tot.append(calculated[key])
	res = sum(tot).doit()
	powe[0] = powe[0].subs(Agabcd,res)
	return powe[0]

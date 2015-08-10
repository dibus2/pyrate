#!/usr/bin/env python 
from RGEsDefinition import *
from ModelsClass import	Model,loggingInfo,loggingDebug,loggingCritical
from RGEsmathModule import DynkinIndex,DynkinCasimir,multiplicity,SimplifyTraces,DeterminOrdering

####################
#Trilinear functions
####################

def set_globalsTri(model) : 
	"""coppy the globals from RGEsModule into the local gloabls()"""
	for key,val in model.glob.items():
		globals()[key] = val 

def CompileTri(Translated,y,comb,model,Weyl,GutNorm):
	"""Compile a given beta function for the Yukawas i.e. divide by the LHS factor"""
	LH = model.TriToCalculate[y][1]
	FinalBeta = (sum([el[0] for el  in Translated[y]]).doit()/LH).expand()
	if Weyl :
		FinalBeta = FinalBeta.subs(kappa,Rational(1,2)) 
	else :
		FinalBeta = FinalBeta.subs(kappa,1) 
	if GutNorm:
            FinalBeta = FinalBeta.subs(model.UGaugeGroups[0][1].g,sqrt(Rational(3,5))*model.UGaugeGroups[0][1].g)
	FinalBeta = FinalBeta.subs(tuple([(el,0) for el in ListAllSymbols['TrilinearTerm']]))
	FinalBeta = DeterminOrdering(model,FinalBeta)
	return FinalBeta

#######
#1-Loop
#######

#power 0 

def CHabc(powe,comb,model):
	"""Calculates the Habc term Eq 66"""
	sc1,sc2,sc3 = comb
	res = model.Expand((('Habc',sc1,sc2,sc3)))
	powe[0] = powe[0].subs(Habc,res)
	return powe[0]

def CL2abc(powe,comb,model):
	"""Calculate the L2abc term Eq 65"""
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {}:
		res = [Rational(1,2)*model.Expand(((_L,a,b,s1,s2),(_h,s1,s2,c)))
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(L2abc,res)
	else :
		powe[0] = powe[0].subs(L2abc,Integer(0))
	return powe[0]

def CLYabc(powe,comb, model):
	"""Calculate the LYabc term. Eq 67"""
	sc1,sc2,sc3 = comb
	res = [model.Expand((('Yab2S',sc,sc),(_h,sc1,sc2,sc3)))
			for sc in [sc1,sc2,sc3]
		]
	res = sum(res)
	powe[0] = powe[0].subs(LYabc,res)
	return powe[0]

#power 2 

def CLSabc(powe,comb,model):
	"""Calculate the term LSabc Eq 68"""
	a,b,c = comb
	res = [model.Expand(((_G,gg1),(_Cs,gg1,s),(_h,a,b,c)))
			for s in [a,b,c]
		]
	res = sum(res)
	powe[0] = powe[0].subs(LSabc,res)
	return powe[0]

######
#2Loop
######

############
#lambda**2*h
############

def CL2habc(powe,comb,model):
	"""Calculates the term sum_L2(i)*habc 69 LINE 1"""
	sc1,sc2,sc3 = comb
	res = [Rational(1,2)*model.Expand((('L2abS',s,s),(_h,sc1,sc2,sc3)))
		for s in [sc1,sc2,sc3]
		]
	res = sum(res)
	powe[0] = powe[0].subs(L2habc,res)
	return powe[0]

def CLbar3abc(powe,comb,model):
	"""Calculates the term Lbar3abc Eq 70"""
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {} :
		res = [Rational(1,2)*(model.Expand(((_L,a,b,s1,s2),(_L,c,s1,s3,s4),(_h,s2,s3,s4)))
			+ model.Expand(((_L,a,s1,s2,s3),(_L,b,s4,s2,s3),(_h,c,s1,s4)))
			)
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(Lbar3abc,res)
	else :
		powe[0] = powe[0].subs(Lbar3abc,Integer(0))

	return powe[0]

#########
#labcd**2
#########

def CLbar2Yabc(powe,comb,model):
	"""Calculates the term CLbar2Yabc Eq 71."""
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {}:
		res = [ Rational(1,2)*model.Expand((('Yab2S',s1,s2),(_L,a,b,s3,s1),(_h,c,s3,s2)))
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(Lbar2Yabc,res)
	else :
		powe[0] = powe[0].subs(Lbar2Yabc,Integer(0))
	return powe[0]

######
#labcd
######

def CHbarlmabc(powe,comb, model):
	"""Calculates the term Hbarlmabc Eq 72"""
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {} : 
		res = [Rational(1,8)*(
			model.Expand(((_L,a,b,s1,s2),('Chain3Ya',s2,c,s1,p1,p2),(_mf,p2,p1)))#Line 1 
			+ model.Expand(((_L,a,b,s1,s2),(_mf,p1,p2),('Chain3Ya',s1,c,s2,p2,p1)))#Line 2 1st term
			+ model.Expand(((_L,a,b,s1,s2),('Chain3Y',s2,c,s1,p1,p2),(_mfa,p2,p1)))#Line 2 2nd term
			+ model.Expand(((_L,a,b,s1,s2),(_mfa,p1,p2),('Chain3Y',s1,c,s2,p2,p1)))
			)
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(Hbarlmabc,res)#norm
	else :
		powe[0] = powe[0].subs(Hbarlmabc,Integer(0))#norm
	return powe[0]

def CH2habc(powe,comb, model):
	"""Calculates the term H2habc Eq 69 3rd line"""
	sc1,sc2,sc3 = comb
	res = [model.Expand((('H2abS',s,s),(_h,sc1,sc2,sc3)))
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(H2habc,res)
	return powe[0]

def CHbar2habc(powe,comb, model):
	"""Calculates the term Hbar2habc Eq 69"""
	sc1,sc2,sc3 = comb
	res =[model.Expand((('Hbar2abS',s,s),(_h,sc1,sc2,sc3)))
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(Hbar2habc,res)
	return powe[0]

def CHbarhabc(powe,comb, model):
	"""Calculates the term Hbarhabc Eq 73"""
	sc1,sc2,sc3 = comb
	res = [Rational(1,4)*(model.Expand(((_h,a,s1,s2),('Chain4Ya',b,s1,c,s2,p1,p1)))
		+ model.Expand(((_h,a,s1,s2),('Chain4Y',s2,b,s1,c,p1,p1)))
		)
		for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
		]
	res = sum(res)
	powe[0] = powe[0].subs(Hbarhabc,res)
	return powe[0]

####
#Yuk
####

def CHYabcHbarYabcH3abc(powe,comb, model):
	"""Calculates the terms HYabc HbarYabc and H3abc to speed up the calculation Eqs 48-49-50"""
	sc1,sc2,sc3 = comb
	tot = []
	Calculated = {}
	for perm in list(permutations([sc1,sc2,sc3],3)):
		key = tuple(flatten(perm))
		if key in Calculated :
			res = Calculated[key]
			tot.append(res)
		else :
			a,b,c = perm
			resHY = (
					model.Expand((('Y2F',p1,p2),(_mfa,p2,p3),('Chain3Y',a,b,c,p3,p1)))#HY Line 1
				+ model.Expand((('Y2F',p1,p2),(_Ya,a,p2,p3),(_mf,p3,p4),(_Ya,b,p4,p5),(_Y,c,p5,p1)))#Line 2
				+ model.Expand((('Y2F',p1,p2),(_Ya,a,p2,p3),(_Y,b,p3,p4),(_mfa,p4,p5),(_Y,c,p5,p1)))#Line 2
				+ model.Expand((('Y2F',p1,p2),('Chain3Ya',a,b,c,p2,p3),(_mf,p3,p1)))#Line 3
				)
			resHYbar = (
					+ 2*Rational(1,2)*(
						model.Expand(((_Y,s1,p1,p2),(_mfa,p2,p3),('Chain4Y',s1,a,b,c,p3,p1)),dotrace=False)#HYbar#Line 1
					+ model.Expand((('Chain4Y',b,c,s1,a,p1,p2),(_Y,s1,p2,p3),(_mfa,p3,p1)),dotrace=False)#Line 2
					+ model.Expand((('Chain3Ya',c,s1,a,p1,p2),(_Y,s1,p2,p3),(_Ya,b,p3,p4),(_mf,p4,p1)),dotrace=False)#Line 3
					+ model.Expand((('Chain4Ya',a,s1,b,c,p1,p2),(_mfa,p2,p3),(_Y,s1,p3,p1)),dotrace=False)#Line 4
					)
					)
			resHYbarad = resHYbar.adjoint() if resHYbar != 0 else 0
			resH3 = (
					+ 2*Rational(1,2)*(#H3abc
					model.Expand(((_mf,p1,p2),('Chain4Ya',a,s1,b,c,p2,p3),(_Ya,s1,p3,p1)))#Line 1
					+ model.Expand(((_mfa,p1,p2),(_Y,s1,p2,p3),('Chain4Ya',b,c,s1,a,p3,p1)))#Line 2
					+ model.Expand((('Chain4Y',c,s1,a,b,p1,p2),(_Y,s1,p2,p3),(_mfa,p3,p1)))#Line 3
					+ model.Expand((('Chain4Y',a,b,s1,c,p1,p2),(_mf,p2,p3),(_Ya,s1,p3,p1)))#Line 4
					)
					)
			resHY = resHY.doit() if resHY != 0 else 0 
			resH3 = resH3.doit() if resH3 != 0 else 0 
			resHYbar = resHYbar.doit() if resHYbar != 0 else 0
			resHYbarad = resHYbar.doit() if resHYbar != 0 else 0 
			res = resHY + SimplifyTraces(Tr(resHYbar+ resHYbarad),model) + resH3
			if res != 0 :
				tot.append(res)
			Calculated[key] = res
	res = sum(tot)
	powe[0] = powe[0].subs(HbarYabc,0).subs(H3abc,0).subs(HYabc,res)
	return powe[0]


#####
#g**2
#####

def CLbar2Sabc(powe,comb,model):
	"""Calculates the term CLbar2Yabc Eq 77."""
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {} :
		res = [Rational(1,2)*model.Expand(((_G,gg1),(_Cs,gg1,s1),(_h,a,s2,s1),(_L,b,c,s2,s1)))
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(Lbar2Sabc,res)
	else :
		powe[0] = powe[0].subs(Lbar2Sabc,Integer(0))

	return powe[0]

def CL2gabc(powe,comb,model):
	"""Calculates the term L2gabc Eq78 """
	sc1,sc2,sc3 = comb
	if model.LbdToCalculate != {} :
		res = [Rational(1,2)*model.Expand(((_G,gg1),(_h,a,s1,s2),(_L,b,c,s3,s4),(_Th,gg1,s1,s3),(_Th,gg1,s2,s4)),MatStruc=['A','A'],dotrace=False)
			for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
			]
		res = sum(res)
		powe[0] = powe[0].subs(L2gabc,res)
	else :
		powe[0] = powe[0].subs(L2gabc,Integer(0))
	return powe[0]


def CY2Fhabc(powe,comb,model):
	"""Calculates the term Y2Fhabc Eq 69 Line 6"""
	sc1,sc2,sc3 = comb
	res = [model.Expand(((_G,gg1),('Y2FabS',gg1,s,s),(_h,sc1,sc2,sc3)))
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(Y2Fhabc,res)
	return powe[0]

def CHSabc(powe,comb, model):
	"""Calculates the term HSabc Eq79 """
	sc1,sc2,sc3 = comb
	res = [model.Expand(((_G,gg1),(_Cs,gg1,s),('Habc',sc1,sc2,sc3)))
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(HSabc,res)
	return powe[0]

def CHFabc(powe,comb, model):
	"""Calculates the term HFabc Eq 80"""
	sc1,sc2,sc3 = comb
	res = [model.Expand(((_G,gg1),(_C,gg1,p1),(_mf,p1,p2),('Chain3Ya',a,b,c,p2,p1)))#Line 1
			+ model.Expand(((_G,gg1),(_mf,p1,p2),(_C,gg1,p2),('Chain3Ya',a,b,c,p2,p1)))#Line 1
			+ model.Expand(((_G,gg1),(_C,gg1,p1),(_Y,a,p1,p2),(_mfa,p2,p3),(_Y,b,p3,p4),(_Ya,c,p4,p1)))#Line 2
			+ model.Expand(((_G,gg1),(_Y,a,p1,p2),(_C,gg1,p2),(_mfa,p2,p3),(_Y,b,p3,p4),(_Ya,c,p4,p1)))#Line 2
			+ model.Expand(((_G,gg1),(_C,gg1,p1),(_Y,a,p1,p2),(_Ya,b,p2,p3),(_mf,p3,p4),(_Ya,c,p4,p1)))#Line 3
			+ model.Expand(((_G,gg1),(_Y,a,p1,p2),(_C,gg1,p2),(_Ya,b,p2,p3),(_mf,p3,p4),(_Ya,c,p4,p1)))#Line 3
			+ model.Expand(((_G,gg1),(_C,gg1,p1),(_Y,a,p1,p2),(_Ya,b,p2,p3),(_Y,c,p3,p4),(_mf,p4,p1)))#Line 4
			+ model.Expand(((_G,gg1),(_Y,a,p1,p2),(_C,gg1,p2),(_Ya,b,p2,p3),(_Y,c,p3,p4),(_mf,p4,p1)))#Line 4
		for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
		]
	res = sum(res)
	powe[0] = powe[0].subs(HFabc,res)
	return powe[0]


#####
#g**4
#####

def Cg4LSabc(powe,comb,model):
	"""Calculate the term LSabc Eq 69 line 7"""
	sc1,sc2,sc3 = comb
	res = [Rational(35,3)*model.Expand(((_G,gg1),(_G,gg1),(_Cg,gg1),(_Cs,gg1,s),(_h,sc1,sc2,sc3)),dotrace=False)
			- Rational(10,3)*kappa*model.Expand(((_G,gg1),(_G,gg1),(_S,gg1),(_Cs,gg1,s),(_h,sc1,sc2,sc3)),dotrace=False)
			- Rational(11,12)*model.Expand(((_G,gg1),(_G,gg1),(_Ss,gg1),(_Cs,gg1,s),(_h,sc1,sc2,sc3)),dotrace=False)
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(g4LSabc,res)
	return powe[0]

def CLSSabc(powe,comb,model):
	"""Calculate the term LSSabc Eq 81"""
	sc1,sc2,sc3 = comb
	res = [model.Expand(((_G,gg1),(_G,gg2),(_Cs,gg1,s),(_Cs,gg2,s),(_h,sc1,sc2,sc3)),dotrace=False)
			for s in [sc1,sc2,sc3]
			]
	res = sum(res)
	powe[0] = powe[0].subs(LSSabc,res)
	return powe[0]

def CAlabc(powe,comb,model):
	"""Calculates the tern Alabc Eq 82"""
	sc1,sc2,sc3 = comb
	res = [(model.Expand(((_G,gg1),(_G,gg2),(_h,a,s1,s2),('Theta4',gg1,gg2,s1,s2,b,c))))
		for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
		]
	res = sum(res)
	powe[0] = powe[0].subs(Alabc,res)
	return powe[0]

def CAbarlabc(powe,comb,model):
	"""Calculates the tern Albarabc Eq 83"""
	sc1,sc2,sc3 = comb
	res = [(model.Expand(((_G,gg1),(_G,gg2),(_h,a,s1,s2),('Theta4',gg1,gg2,b,s1,c,s2))))
		for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
		]
	res = sum(res)
	powe[0] = powe[0].subs(Abarlabc,res)
	return powe[0]

def CBYabcBbarYabc(powe,comb,model):
	"""Calculates the term BYabc and BbarYabc Eqs 84 85"""
	sc1,sc2,sc3 = comb
	Calculated = {}
	tot = []
	for  perm in list(permutations([sc1,sc2,sc3],3)):
		key = tuple(flatten(perm))
		if key in Calculated : 
			res = Calculated[key]
			tot.append(res)
		else :
			a,b,c = perm
			res =(Rational(1,4)*(#BY
		model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_Ta,gg4,p2,p3),(_mf,p3,p4),(_Ya,c,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	+	model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_mf,p1,p2),(_T,gg3,p2,p3),(_T,gg4,p3,p4),(_Ya,c,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	+	model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_Ta,gg4,p2,p3),(_Y,c,p3,p4),(_mfa,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	+	model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Y,c,p1,p2),(_T,gg3,p2,p3),(_T,gg4,p3,p4),(_mfa,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	)
	- 10*Rational(1,4)*(#BYbar
		model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_mf,p2,p3),(_T,gg4,p3,p4),(_Ya,c,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	+	model.Expand(((_G,gg1),(_G,gg2),(_Th,gg1,a,s1),(_Th,gg2,s1,b),(_Ta,gg3,p1,p2),(_Ya,c,p2,p3),(_T,gg4,p3,p4),(_mf,p4,p1)),MatStruc=[['A','B','A','B'],['B','A','A','B']])
	)
	)
			if res != 0 :
				res = res.doit()
				tot.append(res)
			Calculated[key] = res
	res = sum(tot)
	powe[0] = powe[0].subs(BYabc,res).subs(BbarYabc,0)
	return powe[0]

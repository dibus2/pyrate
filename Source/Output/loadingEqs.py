#!/usr/bin/env python 
"""
This contains a function to load a result from a pickle file.
In the Future it is designed to contain a tool box to carry out different operation : 
		- First example, Create a BetaFunction from a well define result stored in a pickle (Do not need to rerun)
"""
import pickle
from sympy import Symbol,transpose
import sys
sys.path.append('../GroupTheory/')
from GroupDefinitions import *
from RGEsmathModule import *
import re as reg

def Load(f):
	#open the file
	Of = open('{}'.format(f),'r')
	try : 
		strres,strsettings,info = pickle.load(Of)
	except :
		exit("problem with the file, either not the correct name or pickle crashed")
	#we first declare all the variables and we also have to declare the classes to make sense of the different expressions
	Replace = []
	for ill,ll in enumerate(strsettings) : 
		#simplify the symbol to get rid of the LateX structure
		lln1 = [''.join(reg.split('{(.*)}',el)) for el in ll if len(reg.split('{(.*)}',el)) == 3 ]
		lln2 = [''.join(reg.split('\\\(.*)',el)) for el in ll if len(reg.split('\\\(.*)',el)) == 3 ]
		lln = lln1 + lln2
		Replace.append([(ll[iel],el) for iel,el in enumerate(lln)])
		#Declare all the symbols
		ll = ll + lln
		for symb in ll : 
			if ill != 0 and ill != 1 : 
				globals()[symb] = Symbol(symb) 
		 	else :
				globals()[symb] = Symbol(symb,commutative = False)
	Replace = sum(Replace,[])
	for rr in Replace :
		strres = [el.replace(*rr) for el in strres]
	res = [eval(el) for el in strres]
	#translate the keys into the symbols
	Nres = {}
	for elem in res :
		for key,vals in elem.items():
			Nres[key] = vals
	return Nres

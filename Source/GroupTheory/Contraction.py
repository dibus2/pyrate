#!/usr/bin/env python
"""
Defines the function needed to do the contraction of different fields"""
from Logging import *
import pudb
try :
	import pickle
	import time
	import re as reg
	import functools,operator
	import copy
	from itertools import permutations,combinations
	import os
        import pudb
except ImportError :
	loggingCritical("Error while loading modules")
try :
	Version = False
	import numpy as np
	vers = np.__version__.split('.')
	if vers[-2] < '5' or (vers[-2] == '5' and vers[-1] < 1):
		Version = True
		raise ImportError
except ImportError :
	if Version :
		loggingCritical("\tnumpy version not compatible please get at least 1.5.1 you have {}.".format('.'.join(vers)),verbose=True)
		exit()
	else :
		loggingCritical("\tError while loading numpy. Check the manual for required modules.",verbose=True)
		exit()
try :
	Version = False
	from sympy import Wild,Symbol,Function,symbols,pi,Rational,zeros,I,sqrt,eye,Matrix,MatrixSymbol,KroneckerDelta,flatten,pprint,IndexedBase,Idx,Integer,Add,Mul,Indexed,Sum,conjugate,adjoint,__version__,Mod
	from sympy.physics.secondquant import evaluate_deltas
	if __version__ != '0.7.2' and __version__ != '0.7.3' :
		Version = True
		raise ImportError
except ImportError :
	if Version :
		loggingCritical("\tsympy version incompatible : {}, please get 0.7.2 (recommended) or 0.7.3 .".format(__version__),verbose=True)
		exit()
	else :
		loggingCritical("\tError while loading sympy. Check the manual for required modules.",verbose=True)
		exit()
localdir = os.path.realpath(os.path.dirname(__file__))
#fdb = open(localdir+'/CGCs.pickle','r')
#fdb = open(localdir+'/CGCstestnorm.pickle','r')
fdb = open(localdir+'/CGCsnew.pickle','r')
db = pickle.load(fdb)
fdb.close()

def GetContractionFactor(dic,Group,CGCs=0):
		"""Seeks the contraction Factor for the contracted particles in dic under the name group.
			Works with any number of fields, at least 2,3,4.
                       F.: Modified February 22 2015 in order to deal with multiple singlets. Actually, this function
                       is called on only one place."""
		#Translate the Qnb into Dynkin labels if given by their dim
		key = []
		for part in dic[Group[0]] :
			key.append(part.Qnb[Group[0]])
		#Get the Factor from the databased
		Match = {2 :'Bilinear',3 :'Trilinear', 4 : 'Quartic'}
		if not(len(key) in Match ):
			loggingCritical("This term is nor a Bilinear, Trilinear nor Quartic : {}\nCheck your model file some term might not be gauge invariant.".format(dic[Group[0]]),verbose=True)
			exit()
		Factor = 0
		if tuple(key) in db[Group[1]._absname][Match[len(key)]] :
                        #F. Check wether it is a list of list
                        if type(db[Group[1]._absname][Match[len(key)]][tuple(key)][0]) == tuple :
                            if CGCs != 0 :
                                loggingCritical("WARNING: `CGCs` specified for an invariant which is unique, ignored!",verbose=True)
    		                Factor = db[Group[1]._absname][Match[len(key)]][tuple(key)]
                            else :
    		                Factor = db[Group[1]._absname][Match[len(key)]][tuple(key)]
                        else :
                            if CGCs == 0 :
                                loggingCritical("ERROR, `CGCs` not specified for the invariant `{}` under {} which as several possible contractions to gauge singlet. Please use the interactive mode (started via the option -idb to check which CGC you want to use.".format(key,Group[1]._absname),verbose=True)
                                exit()
                            else :
                                Factor = db[Group[1]._absname][Match[len(key)]][tuple(key)][CGCs-1]
		else : 
			loggingCritical("The contraction factor for {} {} is not in the db, calculation not implemented yet or the term is not a singlet.".format(key,Group[1]._absname),verbose=True)
			exit()
		#We can create the Function that takes i,j,k,... indices and return the correct contraction factor
		return tuple(Factor)#return the function

def GetDynkinLabel(name,DimR):
	"""Return the Dynkin notation of the irrep of name of dim DimR"""
	#Get the list of possible irrep with this dimension
	try :
		Irreps = db[name]['DimToDynkin'][int(DimR)]
	except KeyError :
		loggingCritical("Error while reading in the databased, {},{} is not in the databased.".format(name,DimR),verbose=True)
		exit()
	if type(Irreps) == list : 
		#exit the Irrep is not uniquely defined, one should use the dynlin label to specified the particles
		loggingCritical("Error the irrep of {}, of dim {}, is not uniquely define, please specify the Dynkin labels in the model file directly.".format(name,DimR),verbose=True)
		exit()
	else :
		return Irreps 


def getdimIrrep(irrep,grp):
	"""returns the dim of the given irrep of the groupname"""
	dim =db[grp._absname]['DynkinToDim'][irrep]
	return dim

def Conjugate(irrep):
	"""return the conjugate irrep """
	if type(irrep) != tuple :
		return -1*irrep
	elif type(irrep) == tuple and len(irrep) == 1 : #SU2 case special treatment 2 and 2bar have the same labels 
		return (irrep[0],True)
	elif type(irrep) == tuple and len(irrep) == 2 and irrep[1] : #SU2 case special treatment 2 and 2bar have the same labels 
		return (irrep[0],)
	elif type(irrep) == tuple and len(irrep) != 1 or len(irrep) == 2 and type(irrep[1]) != bool : 
		return irrep[::-1]
	else :
		loggingCritical("Error, this case is not taken into account in Conjugate of {}.".format(irrep),verbose=True)
		exit()

#############
#Ts functions
#############
class Fabc(Function):

	narg = 7
	is_commutative = True

	@classmethod
	def eval(cls,ff,shape,i1,i2,i3,i4,i5,i6):
		if type(i1) != Symbol and type(i2) != Symbol and type(i3) != Symbol and type(i4) != Symbol and type(i5) != Symbol and type(i6) != Symbol:
			#reconstruct the matrix 
			fabc = np.array(ff).reshape(shape)
			return fabc[i1,i3,i5]*fabc[i2,i4,i6]



class Ts(Function):

	narg = 4
	is_commutative = True

	@classmethod
	def eval(cls,A,mats,f1,f2):
		if type(A) != Symbol and type(f1) != Symbol and type(f2) != Symbol :
			#reconstruct the matrix 
			mat = mats[0][A]
			dim = mats[-1]
			mat = np.array(mat).reshape(dim)
			return mat[f1-1,f2-1]


################
#Factor Function
################

class FF(Function):
	narg = 2
	is_commutative=True

	@classmethod
	def eval(cls, args,contraction):
			"""
			implement generic Factor function
			"""
			if all([el.is_integer for el in args]):
				Components = [el[:-1] for el in contraction]
				Idx = [iel for iel,el in enumerate(Components) if el == args]
				if Idx != [] : 
					assert len(Idx) == 1
					return contraction[Idx[0]][-1]
				else :
					return Integer(0)

class FactorFunction(object):
	"""create a class to store the function
	that calculates the contraction factors"""

	def __init__(self,contraction,Kr=1,indices=[]):
		self.contraction = tuple(contraction)
		#Get the number if indices needed
		self.nbind = len(self.contraction[0])-1#any one would do the trick
		#Get max value for the indices
		self.nmaxind = max([max(el) for el in self.contraction])
		self.Kr = Kr
		self.indices = tuple(indices)
		#replace the deriv indices by the exterior ones
							
	def __repr__(self):
		return "FF"


	def CF(self,indices) : 
		"""creates an instance of the FF class with the corresponding contraction and indices"""
		if all([type(el) == int for el in indices]):
			assert len(indices) == self.nbind
			assert max(indices) <= self.nmaxind
		#For all the indices we have to do the replacement in FF
		f = evaluate_deltas((FF(self.indices,self.contraction)*self.Kr).expand())
		for ie,ind in enumerate(indices) :#matching of the tp indices to the actual ones 
			f = f.subs('deriv{}'.format(ie+1),ind)
		return f

	def update(self,fac):
		"""update the GetContractionFactor method"""
		self.__init__(self.contraction,fac[0],fac[1])

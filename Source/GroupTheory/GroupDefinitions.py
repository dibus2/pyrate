#!/usr/bin/env python
from Contraction import *
""" In this file we define all the symbols and conS*nstants of the groups that are used i.e. SU(N) U(1)"""
m, n,o, p, q ,r ,s ,t,u,v  = map(Wild, ['m','n','o','p', 'q', 'r', 's', 't','u','v' ])

class SUn(object):

	"""Dynamical SU(N) class"""

	def __init__(self,N,name):
		self.name = name
		self._absname = 'SU{}'.format(N)
		self.Algebra = 'A'
		self.N = N
		#dim of the adjoint
		self.d = self.CompDimAdj()
		self.g = symbols('g_{}'.format(self.name))
		self.singlet = '1'
		self.Dynksinglet = tuple((self.N-1)*[0])
		#Repr of the algebra or the fundamental
		self.Matrices = db[self._absname]['Matrices']
		self.HBMat = db[self._absname]['HBmat']
		self.Fond = db[self._absname]['Fond']
		self.Adj = db[self._absname]['Adj']
		#self.fabc = self.GetStructureConstants()
                self.fabc = db[self._absname]['Struc']
		#transformation
		self.fabc = (tuple(flatten(self.fabc)),self.fabc.shape)
		#Convert the  matrices into tuples that can be passed to the symbolic functions
		self.U = False

	def __repr__(self):
		"""Change the representation of SU(N)"""
		return "Instance of SU({}) group".format(self.N)
	
	
	def MatRepr(self,ferm1,ferm2,HB=False):
		"""Get the representation of ferm1,ferm2 """
		#The assertion here is that ferm1 and ferm2 are two different components of the same particle
		if not(HB) :
			assert ferm1 == ferm2
			mat = self.Matrices[ferm1]
#F.sparse			mat = (mat['mat'],mat['shape'])
		else  :
			assert ferm1[1] == ferm2[1] #i.e. same representation
			#read the sequence i.e. PiPi PiSig SigPi 
			mat = self.HBMat[ferm1[1]][(ferm1[0],ferm2[0])]
#F.sparse			mat = (mat['mat'],mat['shape'])
		return mat

	def OuterMatrixProduct(self,M1,M2):
		"""Calculates the outer products of two matrices"""
		p1,q1 = M1.shape
		p2,q2 = M2.shape
		outM = zeros((p1+p2,q1+q2))	
		outM[:p1,:q1] = M1
		outM[p1:,q1:] = M2
		return outM

	def CompDimAdj(self):
		"""return the dim of the Adjoint representation"""
		return self.N**2 -1
		
	def C2(self,irrep):
		""" Casimir for the different irepp of SU(N)
				G stands for the adjoint given by the structure constants."""
		if irrep == self.Dynksinglet:
			return 0
		elif irrep != 'G' :
			return db[self._absname]['Casimir'][irrep]
		elif irrep == 'G' : 
			return db[self._absname]['Casimir'][self.Adj]
		else :
			loggingCritical("Error: no such irrep {}".format(irrep),verbose=True)

	def S2(self,irrep):
		"""Dynkin index for the SU(N) group. For the adjoint we have S2(G)=C2(G) see Langacker p110 for instance"""
		#translate the Dynkin notation
		if irrep == self.Dynksinglet:
			return 0
		elif irrep != 'G' : 
			return db[self._absname]['Dynkin'][irrep]
		elif irrep == 'G' : 
			return db[self._absname]['Dynkin'][self.Adj]
		else : 
			loggingCritical("Error: no such irrep {}".format(irrep),verbose=True)
	
	def sumperso(self,list) : 
		out = list[0]
		for el in list[1:]:
			out += el
		return out

	def GetStructureConstants(self,check=False):
		W = MatrixSymbol('W',self.d,1)
		V = MatrixSymbol('V',self.d,1)
		Vec1 = self.sumperso([W[i,0]*Matrix(self.Matrices[self.Fond]['mat'][i]) for i in range(self.d)])
		Vec2 = self.sumperso([V[i,0]*Matrix(self.Matrices[self.Fond]['mat'][i]) for i in range(self.d)])
		ResTrace =[sum([
			-2*(Matrix(self.Matrices[self.Fond]['mat'][i])[j,k]*KroneckerDelta(l,m) - Matrix(self.Matrices[self.Fond]['mat'][i])[l,m]*KroneckerDelta(k,j))*Vec1[k,l]*Vec2[m,j]
			for j in range(self.N)
			for k in range(self.N)
			for l in range(self.N)
			for m in range(self.N)
			if (Matrix(self.Matrices[self.Fond]['mat'][i])[j,k]*KroneckerDelta(l,m) - Matrix(self.Matrices[self.Fond]['mat'][i])[l,m]*KroneckerDelta(k,j))*Vec1[k,l]*Vec2[m,j] != 0])
			for i in range(self.d)]
		Structures = [[[I*(el.diff(W[i,0])).diff(V[j,0]) for i in range(self.d)] for j in range(self.d)]for el in ResTrace]
		Structures = np.array(Structures).reshape(self.d,self.d,self.d)
		if check :
			Check = [Matrix(self.Matrices[self.Fond]['mat'][i])*Matrix(self.Matrices[self.Fond]['mat'][j]) - Matrix(self.Matrices[self.Fond]['mat'][j])*Matrix(self.Matrices[self.Fond]['mat'][i]) - self.sumperso([I*Structures[i,j,k]*Matrix(self.Matrices[self.Fond]['mat'][k])
					for k in range(self.d)])
					for i in range(self.d)
					for j in range(self.d)
					]
			Check = all([el == zeros((self.N,self.N)) for el in Check])
			if Check:
				return Structures
			else : 
				print "ERROR while determining the structure constants"
		else :
			return Structures

class U1(object):
	"""Defines the U1 group class.
		The Casimir operator of the U1 Group is given by the hypercharge squared over 4. Y**2/4"""

	def __init__(self,name):
		self.d = 1
		self.name = name
		if name != "U1":
			self.g = symbols('g_{}'.format(self.name))
		else :
			self.g = Symbol('g1')
		self._absname = 'U1'
		self.Algebra = ''
		self.singlet = '0'
		self.Dynksinglet = 0
		self.U = True
		self.fabcd = 0#no structure constants

	def __repr__(self):
		"""Change the representation of U1 instances."""
		return "Instance of the U1 group."""
	
	def Structure(self,a,b,c,d):
		return "WARNING U1 factors are abelian i.e. no structure constants"

	def MatRepr(self,f1,f2,HB=False) : 
		"""Representation matrices of U1 group. This is for fermions. For Scalars see HB forms."""
		if not(HB):
			f1Ind,f2Ind = f1[1:],f2[1:]
			Indices = [[f1Ind[i],f2Ind[i]] for i in range(len(f1Ind))]
			return f1[0]*functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Indices],1)
		else : 
			assert f1[1] == f2[1]#hypercharge
			Sequence = (f1[0],f2[0])
			f1Ind,f2Ind = f1[2:],f2[2:]
			#Sanity check
			if not(len(f1Ind) == len(f2Ind)) :
				loggingCritical("Error, inconsistent value detected, contact the author")
				exit()
			#Indices = [[f1Ind[i],f2Ind[i]] for i in range(len(f1Ind))]
			Indices = zip(f1Ind,f2Ind)
			if Sequence == (0,0) : 
				return Rational(1,2)*(f1[1]-f1[1].conjugate())*functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Indices],1)
			elif Sequence == (0,1) : 
				return Rational(1,2)*I*(f1[1]+f1[1].conjugate())*functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Indices],1)
			elif Sequence == (1,0) :
				return Rational(1,2)*-I*(f1[1]+f1[1].conjugate())*functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Indices],1)
			elif Sequence == (1,1) :
				return Rational(1,2)*(f1[1]-f1[1].conjugate())*functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Indices],1)

	def C2(self,hypercharge):
		"""In the Case of U1 the irrep cannot just be F or S, it has to be more precise, for instance it has to match one of the entry of the Particle dictionary defined in the model class."""
		if hypercharge =='G':
			return 0
		else :
			return Rational(hypercharge)**2

	def S2(self,hypercharge):
		"""The dynkin index is equal to the Casimir operator for U1 i.e. the Y**2"""
		if hypercharge == 'G':
			return 0
		else :
			return Rational(hypercharge)**2


#!/usr/bin/env python 
"""We define mathematical functions and tools used for the different calculations."""
from GroupDefinitions import *
def Epsilon(lis):
	"""wrapper for the Eijk function of sympy in order to have it work with a list of arg of arbitrary length"""
	#contruct the string of indices
	lst = ["'{}'".format(l) for l in lis]
	#Construct a string of the code to execute
	st = "Eijk("+','.join(lst)+")"
	return eval(st)


def derivnuc(xpr, field):
	"""does the derivation of a single field"""
	#get the indices and the name of the field by which we take the derivative
	nm,ind = field.args[0],field.args[1:]
	#we have to go through the term and find each term with the same field and for each of them take the derivative
	xprargs = list(xpr.args)#In order to be able to modify it
	Derivated = []
	for xx in xprargs : 
		if type(xx) == Indexed and xx.args[0]==field.args[0] :
			assert len(field.args) == len(xx.args)
			lcxprargs = copy.deepcopy(xprargs)
			temp = [elem for elem in lcxprargs if elem != xx ]
			Kr = [KroneckerDelta(ii,xx.indices[idi]) for idi,ii in enumerate(field.indices) if ii != 0 and xx.indices[idi] != 0 ]
			temp = temp + Kr
			Derivated.append(functools.reduce(operator.mul,temp,1))
	return sum(Derivated)

def derivTensor(xpr,fields):
	""" This is a function that derives an expression containing tensors with indices"""
	if type(xpr) == Add :
		lxpr = xpr.args
	else :
		lxpr = (xpr,)
	if type(fields) == Indexed :
		fields = [fields]
	Final = []#Will be summed
	for eltoderive in lxpr : 
		#Create a Buffer of terms and derive until it is empty
		#The derivation is recursive
		Derived = [eltoderive]
		for ff in fields :
			#This returns a sum of terms all of which have to be derive in turn
			Derived = flatten(Derived)
			TermToDerive = copy.deepcopy(Derived)
			Derived = []
			while TermToDerive != [] :
				toderive = TermToDerive[0]
				TermToDerive = TermToDerive[1:]
				if toderive != 0 :
					res = derivnuc(toderive,ff)
					if type(res) == Add :
						Derived.append(list(res.args))
					else :
						Derived.append([res])
		Final.append(Derived)
	Final = sum(flatten(Final))
	if Final == 0 : 
		return Integer(0)
	else :
		return Final

def ispermutation(llist,perm):
	lclist = [lll for lll in llist]
	for ll in perm :
		if ll in lclist :
			lclist.remove(ll)
		else :
			return False
	if lclist == [] :
		return True

def removeperms(llist):
	lclist = [lll for lll in llist]
	for elem in lclist :
		#get all the permutations of the term
		perms = [el for el in lclist if ispermutation([xx.args[0] for xx in elem],[xx.args[0] for xx in el])]
		for pp in perms[1:]: 
			lclist.remove(pp)
	return lclist


class MatM(Function):

	is_commutative = False
	def __init__(self,arg,*args):
		Function.__init__(self)
		if type(arg[0]) != MatM and type(arg[1]) != MatM : 
			self.arg = arg
			self.indices = args
		elif type(arg[0]) == MatM and type(arg[1]) != MatM :
			self.arg = arg[0].arg+(arg[1],)
			self.indices = args
		elif type(arg[0]) != MatM and type(arg[1]) == MatM : 
			self.arg = (arg[0],)+ arg[1].arg
			self.indices = args

	def update(self):
		return MatM(self.arg,*self.indices)

	def transpose(self):
		self.arg = self.arg[::-1]
		self.arg = tuple([el.transpose() for el in self.arg])
		return MatM(self.arg,*self.indices)

class SP(Function) :
	narg = 2
	is_commutative = True
	def __init__(self,arg1,arg2):
		Function.__init__(self)
		if type(arg1) == MatM and type(arg2) == MatM : 
			self.arg = MatM(arg1.arg),MatM(arg2.arg)
		elif type(arg1) == MatM and type(arg2) != MatM : 
			self.arg = MatM(arg1.arg),arg2
		elif type(arg1) != MatM and type(arg2) == MatM :
			self.arg = arg1,MatM(arg2.arg)
		else :
			self.arg = arg1,arg2
		self.indices = ()

	def update(self) :
		return SP(*self.arg)

class trace(Function) :
	def __init__(self,*args) : 
		Function.__init__(self)
		if len(args) == 1 and type(args[0]) == MatM : 
			self.arg = args[0].arg
		else :
			self.arg = args
		self.is_trace = True
		self.indices = ()
	is_commutative = True
	
	def update(self):
		return trace(*self.arg)
	def conjugate(self):
		return trace(*[el.conjugate() for el in self.arg])
	def adjoint(self):
		return trace(*[el.adjoint() for el in self.arg][::-1])
	def transpose(self):
		return trace(*[el.transpose() for el in self.arg][::-1])

class Tr(Function) :
	def __init__(self,arguments) : 
		Function.__init__(self)
		self.arg = arguments
		self.is_trace = True
		self.indices = ()
	narg = 1
	is_commutative = True

	def conjugate(self):
		return Tr(self.arg.conjugate())
	def adjoint(self):
		return Tr(self.arg.adjoint())
	def transpose(self):
		return Tr(self.arg.transpose())

#We redefine the adjoint so that it can deal with the classes of Yukawas
class adjoint(Function):
	"""
	Conjugate transpose or Hermite conjugation.
	"""
	
	nargs = 1
	
	@classmethod
	def eval(cls, arg):
		obj = arg._eval_adjoint()
		if obj is not None:
			 return obj
		obj = arg._eval_transpose()
		if obj is not None:
			return conjugate(obj)
	
	def _eval_adjoint(self):
		if hasattr(self.args[0],'yuk'):
			return self.args[0].switchindices()
		else :
			return self.args[0]
	
	def _eval_conjugate(self):
		return transpose(self.args[0])
	
	def _eval_derivative(self, x):
		return adjoint(Derivative(self.args[0], x, **{'evaluate': True}))
	
	def _eval_transpose(self):
		return conjugate(self.args[0])
	
	def _latex(self, printer, exp=None, *args):
		arg = printer._print(self.args[0])
		tex = r'%s^{\dag}' % arg
		if exp:
			tex = r'\left(%s\right)^{%s}' % (tex, printer._print(exp))
		return tex
	
	def _pretty(self, printer, *args):
		from sympy.printing.pretty.stringpict import prettyForm
		pform = printer._print(self.args[0], *args)
		if printer._use_unicode:
			pform = pform**prettyForm(u'\u2020')
		else :
			pform = pform**prettyForm('+')
		return pform

	
########################################################################################################
#Other functions 
########################################################################################################

def DynkinIndex(Mod, x, ScalarsOrFermions) : 
	"""Calculates the Dynkin Index for a given group in a given model and for the Scalar or Fermions representation."""
	out = sum([x[1].S2(val.Qnb[x[0]])*multiplicity(x,part,val,Mod) 
			for part,val in Mod.__getattribute__('{}'.format(ScalarsOrFermions)).items()
			if val.Qnb[x[0]] != x[1].Dynksinglet])
	return out

def DynkinCasimir(Mod, x, ScalarsOrFermions):
	"""Calculates the Dynkin times the Casimir for a given group in a given model and for either the scalar or Fermions representation."""
	out =  sum(sum([[
		x[1].g**3*y[1].g**2*y[1].C2(val.Qnb[y[0]])*x[1].S2(val.Qnb[x[0]])*multiplicity(x,part,val,Mod)
		for part,val in Mod.__getattribute__('{}'.format(ScalarsOrFermions)).items()
		if val.Qnb[y[0]]!= y[1].Dynksinglet 
		and val.Qnb[x[0]]!= x[1].Dynksinglet]
		for y in Mod.GaugeGroups],[]))
	return out
	
def multiplicity(x,part,val,Mod) :
	"""Calculate the degeneracy."""
	#We get rid of the - signs for the bar representations and take the quantum numbers of all the groups but the one under which we calculate the RGE and the U1 factors
	Multiplicity = [getdimIrrep(charge,Mod.GetGroupFromName[group])
			for group,charge in val.Qnb.items()
			if Mod.GetGroupFromName[group].name != x[0]
				and not(Mod.GetGroupFromName[group].U)]
	Multiplicity = [el if el != 0 else 1 for el in Multiplicity]#when the dim is zero i.e. singlet the multiplicity should be one
	out = val.Gen*functools.reduce(operator.mul,Multiplicity)
	return out

def SimplifyTraces(expression,Mod) : 
	"""Simplify the traces of the Yukawa, especially gets the minus sign outside"""
	if expression == 0 : 
		return 0
	#Get the list of Yuk from the potential which define the order and add the dagger 
	ListYuk = [Symbol(el) for el in Mod.Potential['Yukawas']]
	#add the list of fermion masses that can also be in the traces
	if 'FermionMasses' in Mod.Potential:
		ListFM = [Symbol(el) for el in Mod.Potential['FermionMasses']]
		ListSymbols = ListYuk + ListFM
	else : 
		ListSymbols = ListYuk
	#copy Matrices with string and ListYuk
	CListYuk = [str(el) for el in ListSymbols]
	def tempfilter1(expression,dic): 
		args = dic[pp].args
		factor = [el for el in args if not(hasattr(el,'arg')) and type(el) != adjoint]
		factor = functools.reduce(operator.mul,factor,1) 
		#Matrices = [el for el in args if type(el)==Symbol or type(el) == adjoint]old version before self.Classes
		Matrices = [el for el in args if hasattr(el,'arg') or (type(el) == adjoint and hasattr(el.args[0],'arg'))]
		#Put the non dagged matrices first in the trace so that we can simplify them
		#calculate the minIndexIn i.e. the first index present in the term according to the ordering of ListYuk
		minIndexIn = min([CListYuk.index(str(el.Symb)) for el in Matrices if hasattr(el,'arg')])
		#Order Matrices according to ListYuk
		while True :
			if type(Matrices[0]) == adjoint :
				out = Matrices[1::]
				out.append(Matrices[0])
				Matrices = out[:]
			elif hasattr(Matrices[0],'Symb') and CListYuk.index(str(Matrices[0].Symb)) == minIndexIn : 
				out = functools.reduce(operator.mul,Matrices)
				break
			else :
				out = Matrices[1::]
				out.append(Matrices[0])
				Matrices = out[:]
		return factor*Tr(out)
	pp = Wild('pp',exclude = (Integer(0),))
	qq = Wild('qq', exclude =(Integer(0),))
	rr = Wild('qq', exclude =(Integer(0),))
	expression = expression.expand()
	temp = expression.replace(lambda xpr : xpr.match(Tr(pp)) != None,
		lambda xpr : tempfilter1(xpr,xpr.match(Tr(pp))))
	temp = temp.subs(Tr(0),0)
	temp = temp.expand()
	return temp


def DeterminOrdering(model,Final):
	"""Determin the ordering of the Structure"""
	if type(Final) != Add : 
		Final = [Final]
	else :
		Final = Final.args
	OrderedTerm = []
	for term in Final :
		args = term.args	
		Out,In,Factor = [],[],[]
		#separtate teh expression in three parts
		for iel,el in enumerate(args):
			if hasattr(el,'is_trace'):
				In.append(el.arg)
			elif hasattr(el,'arg') and not(hasattr(el,'is_trace')) or (type(el) == adjoint and hasattr(el.args[0],'arg')) :
				Out.append(el)
			else :
				Factor.append(el) 
		if Out == In and In == [] : 
			OrderedTerm.append(term)
		else :
			Factor = functools.reduce(operator.mul,Factor,1)
			assert len(In) == 1 or len(In) == 0 
			#Here we have to be carefull
			if len(In) == 1 and hasattr(In[0],'arg') :
				In = In
			elif len(In) == 1 : 
				In = list(In[0].args)
			else :
				In = In
			FOutIndices,FInIndices,ToRemove = [],[],[]
			for iel,el in enumerate(In) :
				el = el.arg if hasattr(el,'arg') else el.args[0].arg
				if model.Particles[str(el[0])].Gen != 1 and model.Particles[str(el[1])].Gen != 1 :
					FInIndices.append((el[0],el[1]))
				elif model.Particles[str(el[0])].Gen != 1 and model.Particles[str(el[1])].Gen == 1 :
					FInIndices.append((el[0],))
				elif model.Particles[str(el[0])].Gen == 1 and model.Particles[str(el[1])].Gen != 1 :
					FInIndices.append((el[1],))
				else :
					#Take the term out of In and multiply the Factor with it since it is a scalar anyway
					tp = In[iel]
					if type(tp) == adjoint :
						Factor = Factor * In[iel].args[0].Symb.conjugate()
					else :
						Factor = Factor * In[iel].Symb
					ToRemove.append(iel)
			#Change the name of the outgoing particles.Add an f for fix index
			In = [el for iel,el in enumerate(In) if not(iel in ToRemove)]
			for idx,el in enumerate(Out) :
				el = el.arg if hasattr(el,'arg') else el.args[0].arg
				if model.Particles[str(el[0])].Gen != 1 and model.Particles[str(el[1])].Gen != 1 :
					if idx == 0 : 
						if len(Out) != 1 :
							FOutIndices.append((Symbol(str(el[0])+'_f'),el[1]))
						else :
							FOutIndices.append((Symbol(str(el[0])+'_f'),Symbol(str(el[1])+'_f')))
					elif idx == len(Out) -1 :
						FOutIndices.append((el[0],Symbol(str(el[1])+'_f')))
					else :
						FOutIndices.append((el[0],el[1]))
				elif model.Particles[str(el[0])].Gen != 1 and model.Particles[str(el[1])].Gen == 1 :
					if idx == 0 : 
						FOutIndices.append((Symbol(str(el[0])+'_f'),))
					else :
						FOutIndices.append((el[0],))
				elif model.Particles[str(el[0])].Gen == 1 and model.Particles[str(el[1])].Gen != 1 :
					if idx == len(Out) - 1 : 
						FOutIndices.append((Symbol(str(el[1])+'_f'),))
					else :	
						FOutIndices.append((el[1],))
				else :
					FOutIndices.append(())
			#Get back only the symbols
			def extractsymbs(lOutorIn):
				cplist = []
				for el in lOutorIn : 
					Transpose = False
					Toadd = 0
					if type(el) != adjoint : 
						if hasattr(el,'arg') : 
							Transpose = False if el.arg == model.DefOrderYuk[str(el.Symb)] else True
						Toadd = el.Symb.transpose() if Transpose else el.Symb
					elif len(model.DimYuk[str(el.args[0].Symb)]) != 0 :
						if hasattr(el.args[0],'arg') : 
							Transpose = False if el.args[0].arg == model.DefOrderYuk[str(el.args[0].Symb)][::-1] else True
						Toadd = el.args[0].Symb.conjugate() if Transpose else el.args[0].Symb.adjoint()#we apply transpose on the adjoint
					else :
						Toadd =  el.args[0].Symb.conjugate()
					cplist.append(Toadd)
				return cplist
			#Old
			#Out = [el.Symb if type(el) != adjoint else (el.args[0].Symb.adjoint() if len(model.DimYuk[str(el.args[0].Symb)]) == 2 else el.args[0].Symb.conjugate()) for el in Out ]
			#In = [el.Symb if type(el) != adjoint else (el.args[0].Symb.adjoint() if len(model.DimYuk[str(el.args[0].Symb)]) == 2 else el.args[0].Symb.conjugate()) for el in In]
			Out = extractsymbs(Out)
			In = extractsymbs(In)
			skip = []
			Comb = [Out,In]
			CombIndices = [FOutIndices,FInIndices]
			Res = []
			for iT,Term in enumerate(Comb) :
				Indices = CombIndices[iT]
				if Term == [] :
					pass 
				else :
					while True : 
						#breaking conditions for in and out terms
						if len(Term) == 1 :
							#Reconstruct the skip terms
							skip = [model.Classes[str(el[0])](el[1]) if str(el[0]) in model.Classes and el[1] != () else 
									(Symbol(str(el[0].args[0]),commutative=True).conjugate() if type(el[0]) == conjugate and el[1] == () else (
										Symbol(str(el[0]),commutative=True) if str(el[0]) in model.Classes and el[1] == () else el[0] )) for el in skip]
							#In case Term is a Yukawa it is necessarily the last one of skip
							if str(Term[0]) in model.Classes :
								if (len(Indices) == 1 and Indices[0] == ()) or Indices == [] :
									Term = Symbol(str(Term[0]),commutative=True)
								elif Indices != [] :
									if len(Indices) == 1 and skip != [] and type(skip[0]) == MatM:
										if Indices[-1][0] == Indices[-1][1] and len(str(Indices[-1][1]).split('_')) == 2 :
											Indices[-1] = Indices[-1][0],Symbol(str(Indices[-1][1]).split('_')[0]+'1_f')
										Term = model.Classes[str(Term[0])](*Indices[-1]).transpose()
									else :
										if Indices[-1][0] == Indices[-1][1] and len(str(Indices[-1][1]).split('_')) == 2 :
											Indices[-1] = Indices[-1][0],Symbol(str(Indices[-1][1]).split('_')[0]+'1_f')
										Term = model.Classes[str(Term[0])](*Indices[-1])
								else :
									loggingCritical("error in determining the product",verbose=True)
									exit()
							#identify traces
							elif len(Indices) == 1 and Indices[0] != () and len(Indices[0]) == 2 :
								if Indices[0][0] == Indices[0][1] and len(str(Indices[0][0]).split('_')) == 1 :
									#it is a trace
									Term = trace(Term[0])
									Term = Term.update()
								elif Indices[0][0] == Indices[0][1] and len(str(Indices[0][0]).split('_')) == 2 :
									#We have to modify the name of the two outgoing particles in this case
									Indices[0] = Indices[0][0],Symbol(str(Indices[0][1]).split('_')[0] + '1_f')
									Term[0].indices = Indices[0]
									Term[0] = Term[0].update()
									Term = Term[0]
								else :
									Term = Term[0]
							elif (len(skip) == 1 and len(Term) == 1
									and hasattr(skip[0],'indices') and hasattr(Term[0],'indices')
									and len(skip[0].indices) == 1 and len(Term[0].indices) == 1 
									and len(str(skip[0].indices[0]).split('_')) != 2
									and skip[0].indices[0] == Term[0].indices[0]):
								#then goddamm it it is a SP 
								Term = SP(skip[0],Term[0])
								Term = Term.update()
								skip = []#empty the skip so that we don't multiply by its content below
							else : 
								Term = Term[0]
							skip0 = functools.reduce(operator.mul,skip,1)
							#At this stage it might be that the skip term is a trace
							if type(Term) == MatM and (len(Term.indices) == 2 and Term.indices[0] == Term.indices[1] and len(str(Term.indices[0]).split('_'))== 1):
								#it is a trace
								Term = trace(*Term.arg)
							#At this stage it might be that the skip term is a trace
							if type(skip0) == MatM and (len(skip0.indices) == 2 and skip0.indices[0] == skip0.indices[1] and len(str(skip0.indices[0]).split('_'))== 1):
								#it is a trace
								skip0 = trace(*skip0.arg)
							Term = skip0*Term
							if type(Term)==list :
								exit("ERROR in determin product")	
							Res.append(Term)
							#empty the skip list 
							skip = []
							break
						#ordering of the terms
						tterm = Term[:2]
						indices = Indices[:2]
						res = analyseterm(tterm,indices,model)
						Term = Term[2:]
						Indices = Indices[2:]
						if res != 'skip':
							Term.insert(0,res)
							if hasattr(res,'indices'):# and res.indices != ():
								Indices.insert(0,res.indices)
						else :
							Term.insert(0,tterm[1])
							if len(indices) != 1 :#if there is only one index left it should not be stored in skip
								skip.append((tterm[0],indices[0]))
								Indices.insert(0,indices[1])
							else :
								skip.append((tterm[0],()))
								Indices.insert(0,indices[0])
			OrderedTerm.append(Factor*functools.reduce(operator.mul,Res,1))
	res = sum(OrderedTerm)
	#replace the symbols
	return res
def analyseterm(term,structure,model) : 
	#takes two Yukawas and return the term with correct ordering 
	#if we are in trace then all the indices are summed while if we are out of trace only the exterior are not
	#When there is a scalar product encounter it can happen that only one index is left
	if len(structure) == 0 : #we are dealing with scalars just return the multiplication
		if type(term[0]) == conjugate and str(list(term[0].atoms())[0]) in model.Classes :
			#term[0] = model.Classes[str(list(term[0].atoms())[0])](list(term[0].atoms())[0]).conjugate()
			term[0] = Symbol(str(list(term[0].atoms())[0]),commutative=True)
		else :
			#term[0] = model.Classes[str(term[0])](term[0]) if str(term[0]) in model.Classes else term[0]
			term[0] = Symbol(str(term[0]),commutative=True)#just to be on the safe side
		if type(term[1]) == conjugate and str(list(term[1].atoms())[0]) in model.Classes :
			term[1] = Symbol(str(list(term[1].atoms())[0]),commutative=True).conjugate()
			#term[1] = model.Classes[str(list(term[1].atoms())[0])](list(term[0].atoms())[0]).conjugate()
		else :
			#term[1] = model.Classes[str(term[1])](term[1]) if str(term[1]) in model.Classes else term[1]
			term[1] = Symbol(str(term[1]),commutative=True)
		return term[0]*term[1]

	if len(structure) == 1 :
		return 'skip'

	if len(structure) == 2 and len(structure[0]) == 0 :
		return 'skip'
	if len(structure) == 2 and len(structure[1]) == 0 : 
		return 'skip'

	###
	#SP and dyadic product
	###

	if len(structure[0]) == 1 and len(structure[0]) == len(structure[1]):
		#Remove the transpose which is useless in this case
		if str(type(term[0])) == 'Dagger' and not(str(type(term[1]))== 'Dagger'):
			return SP(term[0].transpose(),term[1]).update()#Transpose to get rid of the dagger
		elif not(str(type(term[0])) == 'Dagger') and str(type(term[1])) == 'Dagger' : 
			return MatM((term[0],term[1]),structure[0][0],structure[1][0]).update()
		elif type(term[0]) == conjugate and str(type(term[1])) == 'transpose' : 
			return MatM((term[0],term[1]),structure[0][0],structure[1][0]).update()
		elif str(type(term[0])) == 'transpose' and type(term[1]) == conjugate :
			return SP(term[0].transpose(),term[1]).update()
		elif type(term[0]) == MatM and type(term[1]) != MatM :
			if str(type(term[1])) !='Dagger' and str(type(term[1])) != 'transpose'  : 
				return MatM((term[0],term[1].transpose()),structure[0][0],structure[1][0]).update()
			else :
				return MatM((term[0],term[1]),structure[0][0],structure[1][0]).update()
		else :
			return 'skip'


	###########
	#Vector*Mat
	###########
	elif len(structure[0]) == 1 and len(structure[1]) == 2 : #A[s_]B[s_,p_] or A[s_]B[p_,s_]
		if structure[0][0] == structure[1][1] :
			if type(term[0]) == 'Dagger' or str(type(term[0])) == 'transpose':
				out = MatM((term[1],term[0].transpose()),structure[1][0])
			else :
				out = MatM((term[1],term[0]),structure[1][0])	
			#out = MatM((term[1],term[0]),structure[1][0])
		elif structure[0][0] == structure[1][0] :#transpose
			if type(term[0]) == 'Dagger' or str(type(term[0])) == 'transpose':
				out = MatM((term[1].transpose(),term[0].transpose()),structure[1][1])
			else :
				out = MatM((term[1].transpose(),term[0]),structure[1][1])
				#		out = MatM((term[1].transpose(),term[0]),structure[1][1])
		else :
			loggingCritical("Error this is not implemented, contact the authors",verbose=True)
		return out.update()

	###########
	#Mat*Vector
	###########
	elif len(structure[0]) == 2 and len(structure[1]) == 1 : #A[s_,p_]B[s_] or A[s_]B[p_,s_]
		if structure[0][1] == structure[1][0] :#no need to transpose
			out = MatM((term[0],term[1]),structure[0][0])
		elif structure[0][0] == structure[1][0] : #transpose
			out = MatM((term[0].transpose(),term[1]),structure[0][1])
		else :
			loggingCritical("Error this is not implemented, contact the authors",verbose=True)
		return out.update()

	########
	#Mat*Mat
	########
	elif len(structure[0]) == 2 and len(structure[0]) == len(structure[1]) :
		if structure[0][1] == structure[1][0] : 
			out = MatM((term[0],term[1]),structure[0][0],structure[1][1])
		elif structure[0][0] == structure[1][0] : 
			out = MatM((term[0].transpose(),term[1]),structure[0][1],structure[1][1])
		elif structure[0][1] == structure[1][1] :
			out = MatM((term[0],term[1].transpose()),structure[0][0],structure[1][0])
		elif  structure[0][0] == structure[1][1] : 
			out = MatM((term[1],term[0]),structure[1][0],structure[0][1])
		else :
			loggingCritical("Error this is not implemented, contact the authors",verbose=True)
		return out.update()#collapse imbricated statements of MatM
	
	else :
		loggingCritical("Error this is not implemented, contact the authors",verbose=True)
		exit()



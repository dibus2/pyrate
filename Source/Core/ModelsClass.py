import sys
sys.path.append('./Source/GroupTheory')
from GroupDefinitions import *
from RGEsmathModule import *
from Particles import *
import re
try :
	import itertools as itr
	from types import MethodType
except ImportError :
	loggingCritical("Error while loading modules",verbose=True)

class Model(object) :
	"""Takes the Settings read from the YAML file and eventually filtered in RGEs.py and construct the class Model out of it incorporating the groups"""
	def __init__(self,Settings):
		############
		#Definitions
		############
		self._Name = Settings['Name']
		self._Author = Settings['Author']
		self._Date = Settings['Date']
		self.saveSettings = copy.deepcopy(Settings)
		self.GaugeGroups = []
		self.NonUGaugeGroups = []
		self.GaugeGroupsName = []
		self.InvariantResults = {}#for storing the results of the invariants calculated at one loop for the two loop case
		self.Databased = {} #for the results of Y,Lambda, mf ...
		self.DimYuk = {}#store the dimension of the Yukawas
		self.ScalarGroups = {}
		self.IndGaugeGroups = []
		self.GetGroupFromName = {}
		self.getHiggsFromScalar = {}
		self.Particles = {}
		self.Scalars = {}
		self.ScalarsIndex = {}
		self.Fermions = {}
		self.CplxScalars = {}
		self.glob = {}
		self.YukToCalculate = {}
		self.LbdToCalculate = {}
		self.ScMToCalculate = {}
		self.TriToCalculate = {}
		self.FMToCalculate = {}
		self.MixingTerm = {}
		self.ListLbd = []
		self.ListScM = []
		self.ListFM = []
		self.ListTri = []
		self.ToOnly = {}
		self.Combination = {}
		#Validate the entries
		loggingInfo("Validating the model...",verbose=True)
		self.ValidateSettings(Settings)
		loggingInfo("\t\t...done",verbose=True)
		#process the input
		loggingInfo("Processing the model...",verbose=True)
		self.process_inputs(Settings)
		loggingInfo("\t\t...done",verbose=True)
		#define the global symbols
		self.declarations()
		self.ListYukawa = []
		self.savesettingsforoutput = []
		#############################################################################
		#DEF of the class to store the Yuk this is dynamic and hence needs to be here
		#############################################################################
		def constructor(self,*args):
			Function.__init__(self)
			if len(args) == 3 :
				self.Symb = args[0]
				self.arg = args[1:]
			else :
				self.arg = args
			self.yuk = True

		def switchindices(self):
			if hasattr(self,'Symb'): 
				return self.__class__(self.Symb,*self.arg[::-1])
			else :
				return self.__class__(*self.arg[::-1])

		def aadjoint(self) : 
			return adjoint(self.switchindices())
		#############################################################################

		#Extract the Yukawas
		self.Classes = {}
		if 'Yukawas' in self.Potential :
			self.Potential['Yukawas'],self.YukToCalculate,self.ListYukawa = self.ExtractTerm('Yukawas')
		#extract the Quartic terms
		if 'QuarticTerms' in self.Potential: 
			self.Potential['QuarticTerms'],self.LbdToCalculate,self.ListLbd = self.ExtractTerm('QuarticTerms')
		#extract the scalar masses terms
		if 'ScalarMasses' in self.Potential :
			self.Potential['ScalarMasses'],self.ScMToCalculate,self.ListScM = self.ExtractTerm('ScalarMasses')
		###Extract the Fermion mass terms
		if 'FermionMasses' in self.Potential : 
			self.Potential['FermionMasses'],self.FMToCalculate,self.ListFM = self.ExtractTerm('FermionMasses')
		if 'TrilinearTerms' in self.Potential : 
			self.Potential['TrilinearTerms'],self.TriToCalculate,self.ListTri = self.ExtractTerm('TrilinearTerms')
		#Check that there is no symbols that look like M or tM or atM or MatM because it interfers with the output
		for elem in self.ListFM + self.ListYukawa : 
			if elem in ['MatM','atM','tM','M','\lambda','lambda'] : 
				loggingCritical('Change the name of your coupling {}, it interferes with MatM'.format(elem),verbose=True)			
				exit()
               
		# Create the classes for the symbols Yuk and FM
		for listy in self.ListYukawa + self.ListFM: 
			#Create the classes dynamically for each yukawa
			self.Classes[listy] = type(listy,(Function,),{
				'narg':2,'is_commutative': False,'__init__':constructor})
			#Add the above declared method as unbound method of the new class created
			self.Classes[listy].switchindices = MethodType(switchindices, None, self.Classes[listy])
			self.Classes[listy].adjoint = MethodType(aadjoint, None, self.Classes[listy])
		#This is to define wether the Yu is (uRbar, Q) or (Q,uR) used for the vector yukawas
		self.DefOrderYuk = {}
		for key,val in self.YukToCalculate.items():
			self.DefOrderYuk[key] = tuple([el[0] for el in val[-1][0:2]])
		for key,val in self.FMToCalculate.items():
			self.DefOrderYuk[key] = tuple([el[0] for el in val[-1]])
		#Translation dictionary to determin the particle content of the chain in the Expand function
		self.translation = {_Y:['s','f','f'],_Ya:['s','f','f'],_C:['g','f'],_Cs:['g','s'],_Cg:['g'],_S:['g'],_Ss:['g'],_G:['g'],
												_mfa:['f','f'],_mf:['f','f'],_ms:['s','s'],_h:['s','s','s'],_L:['s','s','s','s'],_T:['g','f','f'],_Ta:['g','f','f'],_Th:['g','s','s'],'Y2F':['f','f'],'Y2Fa':['f','f'],'Yab2S':['s','s'],'Chain2Y':['s','s','f','f'],'Chain2Ya':['s','s','f','f'],'Chain3Y':['s','s','s','f','f'],'Chain3Ya':['s','s','s','f','f'],'Chain4Y':['s','s','s','s','f','f'],'Chain4Ya':['s','s','s','s','f','f'],'Chain5Y':['s','s','s','s','s','f','f'],'Chain5Ya':['s','s','s','s','s','f','f'],'Hbar2abS':['s','s'],'H2abS':['s','s'],'Y2FabS':['g','s','s'],'Habcd':['s','s','s','s'],'Theta4':['g','g','s','s','s','s'],'Theta2':['g','s','s','s','s'],'L2abS':['s','s'],'ffabc':['g'],'Hab':['s','s'],'Habc':['s','s','s']}
		
	def __repr__(self):
		"""Change the representation of Model instances"""
		return "This is an instance of the {} model, written by {} on {},\n\t\t Gauge group : {},\n\t\t Contains the different particles :\n\t\tFermions: {}\n\t\tScalars: {}\nPotential :{}\n".format(self._Name,self._Author,self._date,self.GaugeGroups, self.Fermions,self.Scalars,self.Potential)

	def process_inputs(self,Settings):
		"""copy the different inputs from Settings and create the instances of particles"""
		self.GaugeGroupsSettings = Settings['Groups']
		self.NbGroupFactors = len(Settings['Groups'])
		self.get_different_gaugegroups()
		for key,value in Settings.items():
			if key == 'Name':
				self._name = value
			elif key == 'Author':
				self._Author = value
			elif key == 'Date' : 
				self._date = value
			elif key == 'Fermions':
				self.Fermions = value
				#Create the particle and store it in Fermions
				for part,val in value.items():
					self.Fermions[part] = particle(part,val,self.GaugeGroups)
			elif key == 'RealScalars':	
				#Copy the particles in the class
				for part,qnb in value.items():
					Qnb = {'Gen': 1, 'Qnb':qnb}
					self.Scalars[part] = particle(part,Qnb,self.GaugeGroups)
			elif key == 'Potential':
				self.Potential = value
		#Now that the Real Scalars have been created we can create the Cplx one associated
		if 'CplxScalars' in Settings : 
			for part,setts in Settings['CplxScalars'].items():
				setts['Norm'] = self.translatenorm(setts['Norm']) 
				self.CplxScalars[part] = higgsField(part,setts,self.GaugeGroups)
				#declare the real degress of freedom associated to the complex ones
				if not(str(self.CplxScalars[part].RealPart) in self.Scalars) :
					self.Scalars[str(self.CplxScalars[part].RealPart)] = particle(self.CplxScalars[part].RealPart,{'Gen': 1, 'Qnb': self.CplxScalars[part].Qnb},self.GaugeGroups,FromCplx=True)
				if not(str(self.CplxScalars[part].CplxPart) in self.Scalars) :
					self.Scalars[str(self.CplxScalars[part].CplxPart)] = particle(self.CplxScalars[part].CplxPart,{'Gen': 1, 'Qnb': self.CplxScalars[part].Qnb},self.GaugeGroups,FromCplx=True)
		self.Particles.update(self.Fermions)
		self.Particles.update(self.Scalars)
		self.Particles.update(self.CplxScalars)
		#store the names of all the particles
		self.FermionsN = [el._name for el in self.Fermions.values()]
		self.ScalarsN = [el._name for el in self.Scalars.values()]
		for key,val in self.CplxScalars.items():
			for ir,real in enumerate(val.RealComponents):
				if real.args != () :
					kkey = [el for el in real.args if el in self.ScalarsN]
					if kkey == []:
						loggingCritical("error one of the real component is not define in the RealScalars, {}".format(key),verbose=True)
						exit()
					else :
						kkey = kkey[0]
				else :
					kkey = real
				if kkey in self.getHiggsFromScalar : 
					self.getHiggsFromScalar[kkey].append(key) 
				else :
						self.getHiggsFromScalar[kkey] = [key]
						self.ScalarsIndex[kkey] = ir
		#Add the scalars to the getHiggsFromScalar dictionary 
		for key,val in self.Scalars.items() :
			if not(val._name in self.getHiggsFromScalar) : 
				#add the scalar which is not part of any Complex higgs
				self.getHiggsFromScalar[val._name] = [str(val._name)]
				self.ScalarsIndex[val._name] = 0#it is on its own so we put it to zero

	def get_different_gaugegroups(self):
		"""Create the different gauge groups"""
		for GroupName, group in self.GaugeGroupsSettings.items():
			#Identify the group and keep the information of weither or not it is a U1 factor this is a security in case people don't call U1 groups U...
			if group == 'U1' :
				self.GaugeGroups.append([GroupName,U1(GroupName),True])
			elif group.split('U')[0] == 'S':#it a SUn group factor
				self.GaugeGroups.append([GroupName,SUn(int(group.split('U')[1]),'{}'.format(GroupName)),False])
				#add the dim so that we know over which indices to sum later on
				self.IndGaugeGroups.append(self.GaugeGroups[-1][1].N)
			else :
				loggingCritical("ERROR gauge group {} is not known, exiting.".format(group),verbose=True)
				exit()
				#Extract the list of the gauge groups and the corresponding gauge coupling
		self.ListGroups	 = [[self.GaugeGroups[i][1],self.GaugeGroups[i][1].g] for i in range(len(self.GaugeGroups))]
		for name,g,t in self.GaugeGroups :
			if not(t):
				self.NonUGaugeGroups.append([name,g,t])
			self.GaugeGroupsName.append(name)
			self.GetGroupFromName[name] = g
		if self.NonUGaugeGroups == [] : 
			#create a dummy SUN group to avoid problems
			name,g,t = 'SUndum',SUn(2,'SUndum'),False
			self.NonUGaugeGroups.append([name,g,t])
			self.GaugeGroups.append([name,g,t])
			self.GaugeGroupsName.append(name)
			self.GetGroupFromName[name] = g


	def translatenorm(self,normdic):
		"""Translate the norm of a given field"""
		if type(normdic) != str : 
			normdic = Rational(normdic)
			return normdic
		else :
			normdic = normdic.replace('Sqrt','sqrt').replace('i','I')
                        if ('I' in normdic or 'sqrt' in normdic) and not('**' in normdic):
				return eval(normdic)
                        elif '**' in normdic and '/' in normdic :
                            src = re.search('([0-9]{1,2})/([0-9])',normdic)
                            normdic = eval(normdic.replace('{}/{}'.format(src.group(1),src.group(2)),'(Rational({},{}))'.format(src.group(1),src.group(2))))
                            return normdic
			else :
				try :
					return Rational(normdic)
				except :
					exit("Error while determining the norm unknown format {}".format(normdic))

	def areequalqnbs(self,qnba,qnbb):
		"""check if qnba and qnbb are the same or opposite"""
		if type(qnba) != tuple :
			#try if it can be converted into Rational 
			return Rational(qnba) == Rational(qnbb) or Rational(qnba) == - Rational(qnbb)
		if qnba == qnbb :
			return True
		if qnba == qnbb[::-1] :
			return True
		if len(qnba) == 2 and len(qnbb) == 1 and qnba[0] == qnbb[0]:
				return True
		elif len(qnba) == 1 and len(qnbb) == 2 and qnbb[0] == qnba[0]:
			return True
		else :
			return False
	
	def ExtractTerm(self,term):
		"""Extract the Yukawas terms, i.e. expand the terms in component form"""
		###########################
		#Define auxiliary functions 
		###########################
		##
		#Try the different permutations
		##

		def ApplyCGC(Term,combination):
			if type(combination[0]) == tuple :
				combout = [list(ell) for ell in combination]
			else :
				combout = [list(combination)]#in the case there is only one gauge group in the whole contractio it s a tuple
			TempGroup = []
			#The only thing one has to do is to construct the subs and apply them to the factor in Term
			Sub = sum([[(Term[1][iel][iem][0],em) for iem,em in enumerate(el)] for iel,el in enumerate(combout)],[])
			res = Term[-1].subs(Sub)
			return res

		def GetCombInd(CombInd,STOP,Term,counter,Full=False):
			All = False
			tempF = []
			while True :
				try :
					tempcomb = CombInd.next()
					#Construct the permutted tempcomb for each term in Collect
					res = ApplyCGC(Term,tempcomb)
					if res != 0 :
						tempF = [tempcomb,res]
						STOP = True
						break
				except StopIteration :
					counter += 1
					STOP = False
					All = True
					break
			if not(Full) :
				return tempF,STOP,counter
			elif Full and All :
				return tempF,STOP,counter,True
			else :
				return tempF,STOP,counter,False

		def ConstructReturn(el,tempF,Term,Return,ReturnToCalc):
			"""Construct the Return list with the term to calculate """
			LHfactor = tempF[-1]#Left Hand side overall factor by which to divide
			#we need to get the value of each particle under the gauge groups
			particles = [self.Particles[str(ell.args[0])] for ell in Term[0]]
			CpContPart = {}#declare a copy of ContractedParticles
			#update teh ContractedParticles dictionary
			for name,g,U in self.NonUGaugeGroups:
				self.GetContractedParticles(ContractedParticles,name,g,particles)
			outparticles = [[particles[i]._name] for i in range(len(particles))]
			for key,vals in ContractedParticles.items():
				CpContPart[key] = [elval._name for elval in vals]
			NotCharge = [False]
			inm = 0
			for name,g,U in self.NonUGaugeGroups:
				iff = 0
				if all(NotCharge) :
					inm -=1#decrease the group because the previous one was not in tempF
				NotCharge = []
				for ff in outparticles : #NOT SURE THE ORDER IS CONSERVED
					if ff[0] in CpContPart[name]:
						ff.append(tempF[0][inm][iff])
						iff+=1
						NotCharge.append(False)
					else :
						ff.append(0)
						NotCharge.append(True)
				inm+=1#count the group
			return LHfactor,outparticles

		########################
		#End Auxiliary functions
		########################
		translation = {'Yukawas':3,'QuarticTerms':4,'ScalarMasses':2, 'TrilinearTerms':3,'FermionMasses':2}
		################################
		##Test to disentengle the terms
		Skipterms = 0
		################################
		Return = {}
		ReturnToCalc = {}
		SkippedTerms = []
		ListTerm = []
		#Calculate if certain terms mix
		self.MixingTerm[term] = {}
		self.AreThereMixingTerms(term)
		#Need to set the vector under column vectors
		if term == 'Yukawas' or term == 'FermionMasses':
			for iterm,(termkey,termval) in enumerate(self.Potential[term].items()):
				#collect fermions
				if type(termval['Fields'][0]) == list :
					ferms = [[subll for subll in ll if subll in self.Fermions] for ll in termval['Fields']]
					scpos = [[il for il,subll in enumerate(ll) if not(subll in ferms[ill])]for ill,ll in enumerate(termval['Fields'])]
					outfields = [['']*3 for ii in range(len(ferms))] if term =='Yukawas' else [['']*2 for ii in range(len(ferms))]
					if scpos != [] : 
						for idii,ii in enumerate(scpos):
							outfields[idii][-1] = termval['Fields'][idii][ii[0]]
							outfields[idii][0] = ferms[idii][0] if self.Fermions[ferms[idii][0]].Gen != 1 else ferms[idii][1]
							outfields[idii][1] = ferms[idii][[il for il in range(2) if il != ferms[idii].index(outfields[idii][0])][0]]
				else :
					ferms = [ll for ll in termval['Fields'] if ll in self.Fermions]
					scpos = [il for il,ll in enumerate(termval['Fields']) if not(ll in ferms)]
					outfields = ['']*3 if term == 'Yukawas' else ['']*2
					if scpos != [] :
						outfields[-1] = termval['Fields'][scpos[0]]
					outfields[0] = ferms[0] if self.Fermions[ferms[0]].Gen != 1 else ferms[1]
					outfields[1] = ferms[[il for il in range(2) if il != ferms.index(outfields[0])][0]]
					self.Potential[term][termkey]['Fields'] = [ll for ll in outfields]
		for iii,(el,val) in enumerate(self.Potential[term].items()) :
			#In case of mixing we need to have all the terms extract before selecting the right combination
			SKIP = False
			ListTerm.append(el)
			#determin the dimension of the Yukawas two steps : 1) find fermions 2) read their Gen dimension
			if type(val['Fields'][0]) == list :
				if term == 'Yukawas' or term == 'ScalarMasses' or term == 'FermionMasses':
					self.DimYuk[el] =[self.Particles[elem].Gen for elem in val['Fields'][0] if elem in self.Fermions and self.Particles[elem].Gen != 1] 
				#get the quantum numbers
				Ppi = [[self.Particles[subvalfields] for subvalfields in valfields] for valfields in val['Fields']]
			else :
				if term == 'Yukawas' or term == 'ScalarMasses' or term == 'FermionMasses':
					self.DimYuk[el] = [self.Particles[elem].Gen for elem in val['Fields'] if elem in self.Fermions and self.Particles[elem].Gen != 1]
				#get the quantum numbers
				Ppi = [[self.Particles[valfields] for valfields in val['Fields']]]
			#2. Expand the higgs into reals and update the list 
			#TODO CONTRACTION IF HAVE TO FIND A WAY OF CREATING a NEW INSTANCE FOR EACH CONTRACTION I think here would be good
			ExpandedTerm = self.ExpandPotential(Ppi,translation[term],[term,el])
			FinalTerm = []
			TempFinalTerm = {}
			FullSinglet = {}
			for ill,ll in enumerate(ExpandedTerm) : 
				ContractedParticles = {}
				FullSinglet[ill] = False
				for subll in ll : 
					Factor = [] 
					for g in self.NonUGaugeGroups:
						self.GetContractedParticles(ContractedParticles,g[0],g[1],Ppi[ill])
						#Get the contraction factor
						if ContractedParticles[g[0]] != [] : 
                                                        #F. These are modifications for the multiple contraction case
                                                        tpCGCs = 0
                                                        if ('CGCs' in val) and (g[0] in val['CGCs']) :
                                                            if type(val['CGCs'][g[0]]) == int :
                                                                val['CGCs'][g[0]] = [val['CGCs'][g[0]]]
                                                            if not(type(val['CGCs'][g[0]]) == list):
                                                                loggingCritical("Error the `CGCs` option must be an int or a list",verbose=True)
                                                                exit()
                                                            else :
                                                                tpCGCs = val['CGCs'][g[0]][ill]
                                                        Factor.append(GetContractionFactor(ContractedParticles,[g[0],g[1]],CGCs=tpCGCs))
						else :
							Factor.append((0,))#in Order to keep the length of Factor equals to teh length of self.NonUGaugeGroups
                                                #Save Factor for normalizing the CGCs
					if not(all([cc == (0,) for cc in Factor])) :
						#For each group plug the indices in the CGCs
						#The fact that the indices have been generated according to the ordering of self.NonUGaugeGroups and that it s the same for the Factor we now that the indices are in the same order
						indicesff = [0]*len(Factor)
						propindices = [0]*len(Factor)
						for idd,fac in enumerate(Factor) : 
							if fac != (0,) :
							#replace i j is cosmetic
							#NOTE the str in line 342 is a fix because in version 0.7.3 of sympy the method indices returns symbols
								indicesff[idd] = tuple([str(part.indices[idd]).replace('j','i') for part in subll[0] if part.indices[idd] != 0])
								propindices[idd] = [[part.indices[idd],getdimIrrep(self.Particles[str(part.args[0])].Qnb[self.NonUGaugeGroups[idd][0]],self.NonUGaugeGroups[idd][1])] for part in subll[0] if part.indices[idd] != 0]
								Factor[idd] = FF(indicesff[idd],fac)
							else :
								#particles are not charged under this group
								Factor[idd] = 1
							#Add the result of the derivative to the FF function
						#Get rid of the zero in propindices and indicesff if there any
						indicesff = [xx2 for xx2 in indicesff if xx2 != 0]
						propindices = [xx2 for xx2 in propindices if xx2 != 0]
						Factor = evaluate_deltas((functools.reduce(operator.mul,Factor,1)*subll[1]).expand())
						tempkey = subll[0]
                                                #check if the key is a permutation of anyother in TempFinalTerm and if so do the permutations and add it
                                                IsaValidKey = 0,0
                                                if len(ExpandedTerm) == 1 or ill == 0 :
                                                    IsaValidKey = (tempkey in TempFinalTerm),0
                                                else :
                                                    #F.February 2015 In case where there is a sum of terms and the indices are not in 1 to 1 correspondance e.g. (1,2,3)->(2,3,4) (this happes when singlets are in. This function does not work -> Enforce the labeling from 1 to 4 ! see expandpotential function
                                                    IsaValidKey = self.isavalidpermutationkey(tempkey,TempFinalTerm,Factor)
						if IsaValidKey[0]:
							TempFinalTerm[IsaValidKey[1]][-1] = TempFinalTerm[IsaValidKey[1]][-1] + IsaValidKey[-1]
						else :
                                                    assert not(tempkey in TempFinalTerm)
                                                    TempFinalTerm[tempkey] = [subll[0],propindices,Factor]
					else : #Full singlet
						tempkey = subll[0]
						#In case where afull singlet term is given by a zero overall coefficient it should be removed
						if subll[1] != 0 :
							if tempkey in TempFinalTerm :
								TempFinalTerm[tempkey][-1] = TempFinalTerm[tempkey][-1] + subll[1]
							else :
								TempFinalTerm[tempkey] = [subll[0],[],subll[1]]
							FullSinglet[ill] = True
			FinalTerm = TempFinalTerm.values()
			#Sanity check the FullSinglet values should all be the same either True or False
			FullSingletcheck = all([lement == FullSinglet.values()[0] for lement in FullSinglet.values()[1:]])
			if not(FullSinglet):
				loggingCritical("Error, value for `FullSingletcheck` inconsistent, contact the author.")
				exit()
			FullSinglet = FullSinglet.values()[0]
			Return[el] = copy.deepcopy(FinalTerm)
			if self.MixingTerm[term][el] != [] : 
				#Create a list containing all the terms 
				loggingInfo("Warning, some terms are mixing with each other, trying to disentangle them ...",verbose=True)
				SKIP = True
				#save the settings of the loop to be able to finish the work afterwards
				SkippedTerms.append([Return[el],el])
			#One thing that is still do be done at this step is to find a Combination of indices such that the coeff is non zero
			if not(SKIP):
				if Factor != [] and not(FullSinglet):
					counter = 0 
					STOP = False
					while not(STOP) :
						#For each contraction we generate an iterator for the possible contraction
						if counter >= len(Return[el]) :
							exit("error, no non zero contraction found for the term {}".format(el))
						CombInd = itr.product(*[itr.product(*[range(1,x[1]+1) for idx,x in enumerate(element)]) for ielem,element in enumerate(Return[el][counter][1])])
						tempF,STOP,counter = GetCombInd(CombInd,STOP,Return[el][counter],counter)
					if tempF == [] :
						loggingCritical("Error no isospin contraction leading to non zero term were found in {}".format(el),verbose=True)
						exit()
					else : 
						LHfactor,outparticles = ConstructReturn(el,tempF,Return[el][counter],Return,ReturnToCalc)
				elif Factor == [] or FullSinglet:
					#full singlet
					counter = 0 
					LHfactor = Return[el][counter][-1]
					outparticles = [[Symbol(str(ell.args[0]))] + len(self.NonUGaugeGroups)*[0]  for ell in Return[el][counter][0]]
					#fill out the dictionary
				ReturnToCalc[el] = Return[el][counter][0],LHfactor,outparticles

		#Deal with the Mixing Terms
		#1 step : For each mixing term calculate the combination and check that the other ones with which it mixes are zero, if non are found exit
                self.ToOnly[term] = []
                self.Combination[term] = {}
                LastRun = False
		for mixingterm in SkippedTerms : 
			if mixingterm[0] != [] : 
				el = mixingterm[1]
				counter = 0 
				STOP1 = False
				Found = False
				try :
					while not(Found) : 
						#We start by finding a valid combination for the first term 
						while not(STOP1) :
							if counter == len(Return[el]) :
					                    loggingCritical("WARNING, impossible to disentangle the term: {} calculating a combination".format(el),verbose=True)
                                                            self.ToOnly[term].append(el)
                                                            LastRun = True
                                                            SkippedTerms.insert(0,mixingterm)
                                                            Found = True
                                                            break
							CombInd = itr.product(*[itr.product(*[range(1,x[1]+1) for idx,x in enumerate(element)]) for ielem,element in enumerate(Return[el][counter][1])])
							#save cterm
							STOP2 = False
							while not(STOP2) :
								tempF,STOP1,counter,All = GetCombInd(CombInd,STOP1,Return[el][counter],counter,Full=True)#all tells us if all the combinations have been tried
                                                                #Let's initialize it here.
								if All :
									STOP2 = True
								if tempF != [] : 
                                                                        self.Combination[term][el] = []
			                                                Valid = []
									for mixwith in self.MixingTerm[term][el] : 
										#Collect the same terms
										valid = True
										#Find the term which correspond to Return[el] then find the permutation and apply the CGC
										for val in Return[mixwith]:
											#local copy of val
											lcel = [str(elxx.args[0]) for elxx in val[0]]
											for part in [str(elxx.args[0]) for elxx in Return[el][counter][0]] : 
												if str(part) in lcel:
													lcel.remove(str(part)) 
												else :
													break#if one is not in the list it is useless to check the others
											if lcel == [] : 
												#the term is in this term
												TermIdentified = val
												break
										lcparts = [str(elem.args[0]) for elem in Return[el][counter][0]]
										#######
										#F. Bug fix on June 5th 2014: the permutation does not work with multiple group for the scalars
										# Put back zeros for fields not charged under the various gauge groups
										#######
										#Get the structure of the indices 
										tempstructureindices = [elem.indices for elem in Return[el][counter][0]]
										Structureindices = [[elem[i] for elem in tempstructureindices] for i in range(len(self.NonUGaugeGroups))]
										#Get the value of the CGCs
										NewInd = [list(elem) for elem in tempF[0]]
										#substitute the j values by the integers
										newind = []
                                                                                #F. I believe that the list of zeros i.e. non of the particles are charged under a given gauge group are useless and I therefore remove them now (June 6th 2015)
                                                                                Structureindices = [elstruc for elstruc in Structureindices if not(all([elemstruc ==0 for elemstruc in elstruc]))]
										for ielem,elem in enumerate(Structureindices):
											tempnewind = []
											ill = 0
											for ll in elem:
												if ll != 0 :
													tempnewind.append(NewInd[ielem][ill])
												else :
													ill -= 1
													tempnewind.append(0)
												ill+=1
											newind.append(tempnewind)
										NewInd = newind
										lctterm = [str(elem.args[0]) for elem in TermIdentified[0]]
										while True :
											NotAllowed = []
											for ix,x in enumerate(lcparts) : #go through the not ordered list and find the first place where it should be 
												if ix in NotAllowed :
													pass
												else :
													for new,nn in enumerate(lctterm) :
														if new in NotAllowed :
															pass
														else :
															if x == nn and not(new in NotAllowed) :
																if ix == new : 
																	NotAllowed.append(ix)
																	break
																else :
																	if lcparts[ix] != lcparts[new] :
																		lcparts[ix] = lcparts[new]
																		lcparts[new] = x
																		#permut the indices for all the groups
																		for idjj,jj in enumerate(NewInd) :
																			buff = jj[ix]
																			jj[ix] = jj[new]
																			jj[new] = buff
																		NotAllowed.append(ix)
																		NotAllowed.append(new)
																		break
																	else :
																		break
											if lctterm == lcparts:
												break
											else :
												pass
										#Get rid of the zeros again
										NewInd = [[elem for elem in ll if elem != 0] for ll in NewInd]
										#Calculate the permutation with respect to the first list
										res = ApplyCGC(TermIdentified,tuple([tuple(lll) for lll in NewInd]))
										if res != 0 : 
											#it is not a valid combination
                                                                                        if LastRun :
                                                                                            pass
                                                                                        else :
                                                                                            self.Combination[term][el].append((mixwith,res))
											    Valid.append(False)
    										    	    break
										if valid : 
											Valid.append(True)
										else :
											STOP = False
                                                                                        self.Combination[term][el].append((mixwith,res))
											break
                                                                                self.Combination[term][el].append((mixwith,res))
									if all(Valid) and Valid != []:
										Found = True
										STOP2= True
                                                                                LastRun = False
										LHfactor,outparticles = ConstructReturn(el,tempF,Return[el][counter],Return,ReturnToCalc)
										ReturnToCalc[el] = Return[el][counter][0],LHfactor,outparticles		
				except IOError :
					exit("error while disentengling should not get there, contact the authors")
					pass
		#Just before exiting the function one has to devide all the combinations by the overall factor
		for key2,val2 in self.Combination[term].items():
			for iv,vv in enumerate(val2):
				self.Combination[term][key2][iv] = vv[0],vv[1]/ReturnToCalc[key2][1]
		#F. March 9th 2015 
		#One has to extract the contribution comming from the CGCs such that it can be normalized to one. This is important since the normalization of the CGCs since Susyno v3 are different.
		#This has to be done for each one of the terms
		#for ll in ListTerm:
		#		#Get the corresponding term from the potential
		#		tppotentialterm = [el[-1] for iel,el in enumerate(Return[ll]) if el[0] == ReturnToCalc[ll][0]]
		#		assert len(tppotentialterm) == 1
		#		tppotentialterm = tppotentialterm[0]
		#		#split up the contributions and keep only the CGCs i.e. FF terms
		#		coefficient = np.array([functools.reduce(operator.mul,[ell for ell in el.args if type(ell) != FF],1) for el in tppotentialterm.args])
		#		coefficient = coefficient*1/coefficient[0]
		#		tppotentialterm = [functools.reduce(operator.mul,[ell for ell in el.args if type(ell) == FF],1)*norm for el,norm in zip(tppotentialterm.args,coefficient)]
		#		#Create list of substitution
		#		tosubs = tuple(flatten([tuple(zip(ff.indices,gg[1:])) for ff,gg in zip(ReturnToCalc[ll][0],ReturnToCalc[ll][-1])],1))
		#		CGCsnorm = np.array([el.subs(tosubs) for el in tppotentialterm])
		#		#At this stage all the different CGCs should return the same overal factor 
		#		checkCGCsnorm = np.all(abs(CGCsnorm)-abs(CGCsnorm[0])==0)
		#		#TODO The pb is how do we know which overall factor to use ? Plus or minus sign? I think the one of the selected combination
		#		if not(checkCGCsnorm) :
		#				loggingCritical("Error while determining the overall normalization of the CGCs, please contact the author.",verbose=True)
		#		tpaddtoReturnToCalc = list(ReturnToCalc[ll])
		#		tpaddtoReturnToCalc.append(CGCsnorm[0])
		#		print(CGCsnorm,tpaddtoReturnToCalc)
		#		tpaddtoReturnToCalc[1] = tpaddtoReturnToCalc[1]/abs(CGCsnorm[0])
		#		ReturnToCalc[ll] = tuple(tpaddtoReturnToCalc)
		#		for iel,el in enumerate(Return[ll]):
		#				Return[ll][iel] = [el[0],el[1],el[2]*1/abs(CGCsnorm[0])]
		return Return,ReturnToCalc,ListTerm


#################################
#aux function for the constructor
#################################
		
	def GetContractedParticles(self,ContractedParticles,name,g,fields):
		"""fill up the ContractedParticles dict with the particles that are contracted together"""
		#let's find which particles are contracted with each other
		ContractedParticles[name] = []
		for ip,part in enumerate(fields) : 
			if part.Qnb[name] != g.Dynksinglet : 
				ContractedParticles[name].append(fields[ip])

	def ValidateSettings(self,Settings):
		"""Implements the different checks carried out on the input provided by the user"""
		#0.1 Check the Mandatory fields are in the model File and that teh particles are consistant
		if not('Groups') in Settings :
			loggingCritical("error no Gauge Groups specified, exit.",verbose=True)
			exit()
		else :
			Groups = Settings['Groups'].keys()
		#Copy the Fermions and Scalars if present
		if 'Fermions' in Settings :
			Fermions = Settings['Fermions']
		else :
			Fermions = {}
		if 'RealScalars' in Settings :
			Scalars = Settings['RealScalars']
		else :
			Scalars = {}
		if 'CplxScalars' in Settings :
			CplxScalars = Settings['CplxScalars']
		else :
			CplxScalars = {}
		TempParticles = {}
		TempParticles.update(Fermions)
		TempParticles.update(Scalars)
		TempParticles.update(CplxScalars)
		#Check that all the gauge groups are define above
		for part,val in TempParticles.items() :
			if 'Qnb' in val :
				Tags = val['Qnb'].keys()
			else :
				Tags = val.keys()
			Tags = all([el in Groups for el in Tags])
			if not(Tags) : 
				loggingCritical("Error, the particle {} is charged under an unknown gauge group.".format(part),verbose=True)
				exit()
		if not('Potential') in Settings :
			loggingCritical("error no Potential specified, exit.",verbose=True)
			exit()
		for key,term in Settings['Potential'].items() :
			Len = {'QuarticTerms':4,'TrilinearTerms':3,'ScalarMasses':2}
			if key == 'QuarticTerms' or key == 'ScalarMasses' or key == 'TrilinearTerms':
				for qq,xpr in term.items():
					#all the fields should be scalars cplx or real
					if not(type(xpr['Fields'][0]) == list) :
						templist = [xpr['Fields']]
					else :
						templist = [el for el in xpr['Fields']]
					for subxpr in templist :
						Tags = [el in Scalars.keys() or el in CplxScalars.keys() for el in subxpr]
						Tags = len(Tags) == Len[key] and all(Tags)
						if not(Tags) :
							loggingCritical("Error, the {} {} does not contain only scalars or the wrong number".format(key,subxpr),verbose=True)
							exit()
			if key == 'Yukawas' :
				for yy,xpr in term.items():
					#two fermions and a scalar field
					if not(type(xpr['Fields'][0]) == list) :
						templist = [xpr['Fields']]
					else :
						templist = [el for el in xpr['Fields']]
					for subxpr in templist :
						Scalar = [el for el in subxpr if el in Scalars.keys() or el in CplxScalars.keys()]
						Fermion = [el for el in subxpr if el in Fermions.keys()]
						Tags = len(Scalar) == 1 and len(Fermion) == 2
						if not(Tags) : 
							loggingCritical("Error, the yukawa term {} does not contain two fermions and a scalar: {}".format(yy,subxpr),verbose=True)
							exit()
			if key == 'FermionMasses':
				for ff,xpr in term.items():
					#two fermions
					if not(type(xpr['Fields'][0]) ==list) :
						templist = [xpr['Fields']]
					else :
						templist = [el for el in xpr['Fields']]
					for subxpr in templist :
						Tags = len([el for el in subxpr if el in Fermions]) == 2
						if not(Tags):
							loggingCritical("Error, the FermionMasses {} does not contain two fermions : {}".format(ff,subxpr),verbose=True)
							exit()

		for key,term in Settings['Potential'].items() :
                    for kkey,tterm in term.items():
                        if 'CGCs' in tterm : 
                            if type(term[kkey]['Fields'][0]) == list : 
                                if not(all([type(valel) == list for keyel,valel in tterm['CGCs'].items()])) or not(all([len(valel)==len(term[kkey]['Fields']) for keyel,valel in tterm['CGCs'].items()])):
                                    loggingCritical("Error, `CGCs` keys must be a dictionary of list of the same size as the `Field` key.",verbose=True)
                                    exit()

                #            if not(all([el in Settings['Groups'].keys() for el in tterm['CGCs']])):
                #                loggingCritical("Error, `CGCs` keys must be known gauge group: {}".format(tterm['CGCs']),verbose=True)
                #                exit()
                #            if any([type(valgr) != int  for kgroup,valgr in tterm['CGCs'].items()]):
                #                loggingCritical("Error, `CGCs` entries must be int: {}".format(tterm['CGCs']),verbose=True)
                #                exit()


	def AreThereMixingTerms(self,term):
		"""
		Checks if they are mixing terms and if yes returns the mixing terms.
		F. 28/01/2013
		In the case where some terms are sums of terms we are assuming that we have 
		checked already that they are consistent i.e. they are composed of the same real
		degrees of freedom.
		"""
		lcPotentialFields = [self.Potential[term].items()[ii][1]['Fields'][0]
				if type(self.Potential[term].items()[ii][1]['Fields'][0]) == list
				else self.Potential[term].items()[ii][1]['Fields']
				for ii in range(len(self.Potential[term]))]
		Fields = [[[self.Particles[f].RealPart,self.Particles[f].CplxPart]
			if self.Particles[f].Cplx 
			else [self.Particles[f]._name] 
			for iff,f in enumerate(lcPotentialFields[ii])] for ii in range(len(lcPotentialFields))
			]
		for iff,ff in enumerate(Fields) : 
			self.MixingTerm[term][self.Potential[term].keys()[iff]] = []
			for igg,gg in enumerate(Fields) :
				if iff != igg :
					lcFields = copy.deepcopy(gg)
					for ffelem in ff :
						if ffelem in lcFields : 
							lcFields.remove(ffelem)
						else :
							break
						if lcFields == [] :
							self.MixingTerm[term][self.Potential[term].keys()[iff]].append(self.Potential[term].keys()[igg])
		return self.MixingTerm 
	

	def ExpandPotential(self,fields,nbfields,term):
		"""Expand the Potential	"""
		SumFOut= []
		#Replace the Complex higgs and create IndexedBase object
		#Determin which configurations have to be calculated, i.e., which derivatives we will take
		if not(type(self.Potential[term[0]][term[1]]['Norm'])) ==list :
			lcNorm = [self.Potential[term[0]][term[1]]['Norm']]
		else :
			lcNorm = self.Potential[term[0]][term[1]]['Norm']
		for idelem,elemfields in enumerate(fields):
			Out,ToDerive,FOut = [],[],[]
                        #F. February 25th There is actually a pretty big flaw here. I need different index count on group basis to skip singlets!!
                        IndicesCounters,dummyff=[0]*len(self.NonUGaugeGroups),0
			for f in elemfields:
				Singlet = True if all([f.Qnb[g[0]] == g[1].Dynksinglet for g in self.NonUGaugeGroups]) else False
				if f.Cplx and not(Singlet): 
					Out.append([IndexedBase(f.RealPart)[['i{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]] + f.Coeff*IndexedBase(f.CplxPart)[['i{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]],f.norm])
					ToDerive.append([IndexedBase(f.RealPart)[['j{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]],IndexedBase(f.CplxPart)[['j{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]]])
                                        #Increment the proper indices
				elif f.Cplx and Singlet :
					#If self.NonUGaugeGroups is an empty list it crashes F. on the 22.07.14
					Out.append([IndexedBase(f.RealPart)[[Symbol('dumi{}{}'.format(dummyff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]] + f.Coeff*IndexedBase(f.CplxPart)[[Symbol('dumi{}{}'.format(dummyff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]],f.norm])
					ToDerive.append([IndexedBase(f.RealPart)[[Symbol('dumj{}{}'.format(dummyff+1,igg))  for igg,g in enumerate(self.NonUGaugeGroups)]],IndexedBase(f.CplxPart)[[Symbol('dumj{}{}'.format(dummyff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]]])
                                        dummyff+=1
				elif not(f.Cplx) and not(Singlet) :
					Out.append(IndexedBase(f._name)[['i{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]])
					ToDerive.append([IndexedBase(f._name)[['j{}{}'.format(IndicesCounters[igg]+1,igg) if f.Qnb[g[0]] != g[1].Dynksinglet else Integer(0) for igg,g in enumerate(self.NonUGaugeGroups)]]])
				else :
					#If self.NonUGaugeGroups is an empty list it crashes F. on the 22.07.14
					#Out.append(IndexedBase(f._name)[[Symbol('dumi{}{}'.format(iff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]])
					#ToDerive.append([IndexedBase(f._name)[[Symbol('dumj{}{}'.format(iff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]]])
					Out.append(IndexedBase(f._name)[[Symbol('dumi{}{}'.format(dummyff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]])
					ToDerive.append([IndexedBase(f._name)[[Symbol('dumj{}{}'.format(dummyff+1,igg)) for igg,g in enumerate(self.NonUGaugeGroups)]]])
                                        dummyff+=1
                                Toinc = [True if f.Qnb[g[0]] != g[1].Dynksinglet else False for igg,g in enumerate(self.NonUGaugeGroups)]
                                IndicesCounters = [el+1  if Toinc[iel] else el for iel,el in enumerate(IndicesCounters)]
			#Take the cartesian product and remove the permutations from the list
			ToDerive = removeperms(list(itr.product(*ToDerive)))
			##If there're several singlets we are screwed so we need to declare dummy indices for them
			Out = (self.translatenorm(lcNorm[idelem])*functools.reduce(operator.mul,[functools.reduce(operator.mul,el,1)  if type(el) == list else el for el in Out],1)).expand()
#                        FullNormNoCGCs = [ell.args for ell in Out.args]
#                        FullNormNoCGCs = [functools.reduce(operator.mul,[el for el in ell if type(el) != Indexed],1) for ell in FullNormNoCGCs]
			#Now we can sum the terms and derive them all together
			for toderive in ToDerive : 
				FOut.append([toderive,derivTensor(Out,toderive)])
			##At this point we have all the derivatives for this term
			CollectDummy = [[(el,Integer(0)) for el in elem[1].atoms() if len(str(el).split('dum')) != 1] for elem in FOut]
#			SumFOut.append([[[tuple([xx.subs(tuple(CollectDummy[iel])) for xx in el[0]]),el[1].subs(tuple(CollectDummy[iel]))] for iel,el in enumerate(FOut)],FullNormNoCGCs])
			SumFOut.append([[tuple([xx.subs(tuple(CollectDummy[iel])) for xx in el[0]]),el[1].subs(tuple(CollectDummy[iel]))] for iel,el in enumerate(FOut)])
		return SumFOut
        
        
        def isavalidpermutationkey(self, term, dictkeys,Factor): 
            """
            determines if term is a permutation of any key in dictkeys. 
            If so returns the list of replacement to do to obtain the new coefficient
            Indices can be repeated
            """
            lcdictkeys = dictkeys.keys()
            TermIdentified = 0
            Baseterm = [xx.base for xx in term]
            Indterm = [xx.indices for xx in term]
            for ll1 in lcdictkeys :
                lcll1 = [xx.base for xx in ll1]
                for ll2 in Baseterm : 
                    if ll2 in lcll1 : 
                        lcll1.remove(ll2)
                    else :
                        break
                if lcll1 == [] :
                    #There is only one possible solution
                    TermIdentified = ll1
                    break
            if TermIdentified == 0 :
                return False,[],Factor
            else :
                #determin the permutations needed
                NewInd = [xx.indices for xx in ll1]
                lcll1 = [xx.base for xx in ll1]
                while True :
                    NotAllowed = []
                    for ix,x in enumerate(Baseterm) :
                        if ix in NotAllowed :
                            pass
                        else : 
                            for new,nn in enumerate(lcll1):
                                if new in NotAllowed:
                                    pass
                                else :
                                    if x == nn and not(new in NotAllowed):
                                        if ix == new :
                                            NotAllowed.append(ix)
                                            break
                                        else :
                                            if Baseterm[ix] != Baseterm[new] : 
                                                Baseterm[ix] = Baseterm[new]
                                                Baseterm[new] = x
                                                Factor = self.swapindices(Factor,NewInd,ix,new)
                                                buff = NewInd[ix]
                                                NewInd[ix] = NewInd[new]
                                                NewInd[new] = buff
                                                NotAllowed.append(ix)
                                                NotAllowed.append(new)
                                                break
                                            else  :
                                                break
                    if lcll1 == Baseterm : 
                        break
                    else :
                        pass 
                return True,TermIdentified,Factor

        def swapindices(self,Factor,NewInd,ix,new):
            "Swap indices ix,new in Factor"
            for idx, xx in enumerate(NewInd[ix]) :
                if NewInd[new][idx] != 0 and xx != 0 :
                    Factor = Factor.subs(xx,'temp')
                    Factor = Factor.subs(NewInd[new][idx],xx)
                    Factor = Factor.subs('temp',NewInd[new][idx])
            return Factor

##				for tterm in TermIdentified : 
#					lcparts = [str(elem) for elem in parts]
#					NewInd = [elem for elem in indices]
#					lctterm = [str(elem.args[0]) for elem in tterm[1][0]]
#					while True :
#						NotAllowed = []
#						for ix,x in enumerate(lcparts) : #go through the not ordered list and find the first place where it should be 
#							if ix in NotAllowed :
#								pass
#							else :
#								for new,nn in enumerate(lctterm) :
#									if new in NotAllowed :
#										pass
#									else :
#										if x == nn and not(new in NotAllowed) :
#											if ix == new : 
#												NotAllowed.append(ix)
#												break
#											else :
#												if lcparts[ix] != lcparts[new] :
#													lcparts[ix] = lcparts[new]
#													lcparts[new] = x
#													buff = NewInd[ix]
#													NewInd[ix] = NewInd[new]
#													NewInd[new] = buff
#													NotAllowed.append(ix)
#													NotAllowed.append(new)
#													break
#												else :
#													break
#						if lctterm == lcparts:
#							break
#						else :
#							pass
   




##########################
#End constructor functions
##########################

#################
#Expand functions
#################

	def Expand2(self,args,dotrace=True,Layer=0,MatStruc=[]):
		"""Expand the term with valid particles
			Each term is given by ('Y',p1,p2,s1) ...
			MatStruc is used for the Mat representations. It gives the structure of the chain in this case.
			Layer counts the nested calls to Expand needed for the saved versions
			"""
		Layer += 1
		#Check if there is a Trilinear,scalar mass, fermion mass
		argsflt = flatten(args)
		TriIn = _h in argsflt
		FMIn = _mf in argsflt or _mfa in argsflt
		ScMIn = _ms in argsflt
		if TriIn and self.TriToCalculate == {} :
			return 0
		if FMIn and self.FMToCalculate == {} :
			return 0
		if ScMIn and self.ScMToCalculate == {} :
			return 0
		#we start by identifying all the particles that go together
		scalars = self.ScalarsN
		fermions = self.FermionsN
		groups = self.GaugeGroupsName 
		#check if there is only one argument transform the tuple
		if type(args[0]) != tuple and type(args[0]) != list :
			args = (args,)
		if type(args) == list : 
			args = [args]
		if dotrace :
			#we ignore the substructure irrelevant for determining the trace
			argsSimple = [el[0] for el in args]
			argsS,Structure = [],[]
			for el in args : 
				if type(el) == list : 
					for ell in el[0] : 
						argsS.append(ell)
						Structure.append(ell[0])
				else :
					argsS.append(el)
					Structure.append(el[0])
			Traces = self.DeterminTraces2(argsS,Structure)
		else : 
			Traces = (None,None)
		###################
		#end Initialization
		###################
		#determin the structure in terms of blocks e.g. A(p1,p2,p3)B(p3,p4,p5)C(p5,p6,p1) 
		Blocks = []
		temp = []
		for el in args : 
			if type(el) != list : 
				temp.append(el)
			else : 
				#store temp
				if temp!= []:
					Blocks.append(temp) 
					temp = []#reset it to zero
				Blocks.append(el)#store the new block
		if temp != [] :
			Blocks.append(temp)
		#Determin the structure of the chain
		SaveChain = []
		TChain = []
		FullTChain = []
		FixedIndices = []
		BlockMapping = {}
		for ielem,elem in enumerate(Blocks):
			if type(elem[0]) == list : #sum of terms one more layer to go down
				SaveChain.append(sum([list(el[1:]) for el in elem[0]],[]))#the structure is the same in all the elements of one block no need to go through 
				FullTChain.append([[el[0] for el in ii] for ii in elem])
				TChain.append([el[0] for el in elem[0]])
				BlockMapping[ielem] = {}
				#we need to get the mapping of the indices so taht we can coppy back the indices once generated see below
				for ibbl,bblock in enumerate(elem[1:]):
					tpbblock = sum([list(x[1:]) for x in bblock],[])
				#	alreadycopied = []
					BlockMapping[ielem][ibbl] = []
					for ieeb,eeb in enumerate(tpbblock) : 
						for iib,eb in enumerate(SaveChain[-1]):
							if eb == eeb:# and not(iib in alreadycopied) : 
								BlockMapping[ielem][ibbl].append(iib)
								#alreadycopied.append(iib)
								break
			else : #tuples product terms 
				SaveChain.append(sum([list(el[1:]) for el in elem],[]))
				FullTChain.append([el[0] for el in elem])
				TChain.append([el[0] for el in elem])
		#we need to identify the indices that are fixed in the chain. These appear as non global variables
		SumSaveChain = sum(SaveChain,[])
		SumTChain = flatten(TChain)
		FixedIndices = [iel for iel,el in enumerate(SumSaveChain) if type(el) == list or not(str(el) in self.glob)]
		#Link between SumSaveChain and SaveChain i.e. gives how many indices go in each block
		#LinkSumSaveChain = [len(el) - [True if type(ell) == list or not(str(ell) in self.glob) else False for ell in el].count(True) for el in SaveChain] 
		LinkSumSaveChain = [len(el)  for el in SaveChain] 
		#remove the fix indices from the generation
		Chain = [el for iel,el in enumerate(SumSaveChain) if not(iel in FixedIndices)]
		#translate using the dictionary above
		TypeChain = [self.translation[el] for el in SumTChain]
		TypeChain = sum(TypeChain,[])
		TypeChain = [el for iel,el in enumerate(TypeChain) if not(iel in FixedIndices)]
		#do the link between the indices in the chain and the generated indices below
		GeneratePosition = []
		RemovedIndices = []
		PartInd = {}
		ip = 0
		for it,term in enumerate(Chain) :
			if term in PartInd:
				ind1 = PartInd[term]#index of the last identic term
				#found the corresponding value in Generate position
				ind2 = [ii[1] for ii in GeneratePosition if ii[0] == ind1]
				assert len(ind2) ==1
				GeneratePosition.append((it,ind2[0]))
				RemovedIndices.append(it)
			else :
				PartInd[term] = it#create the entry
				GeneratePosition.append((it,ip))
				ip+=1
		#create the lists for the cartesian product
		Lists = [fermions if el == 'f' else (scalars if el == 's' else groups) 
				for iel,el in enumerate(TypeChain) if not(iel in RemovedIndices)]
		#CopyLists = list(Lists)
		CopyLists = copy.deepcopy(Lists)
		#create the iterator i.e. cartesian product
		it = itr.product(*Lists)
		FinalTerm = []#for the final result
		#iteration over the combination
		for icm,combo in enumerate(it) : 
			ExtraConstrains = [] #used to do some matching for the Theta and t matrices
			STOP = False
			GenIndices = [combo[el[1]] for el in GeneratePosition]
			#Indices = 
			Out = []
			iigen = Index(self)#create a generator of indices
			#generate the indices qui vont bien 
			indices = iigen.GenConsInd(combo,Layer)
			#duplicate the indices that are repeated
			indicesFull = [indices[el[1]] for el in GeneratePosition]#position in Generate are given with respect to the Chain list :w
			#SummedIndices = list(indices)
			SummedIndices = copy.deepcopy(indices)
			#add the ones that are fixed in the correct place
			for el in FixedIndices : 
				extindices = iigen.GenExtInd(SumSaveChain[el],Layer)
				if type(SumSaveChain[el]) == str :#group factors are not given under a list but just the name
					partname = SumSaveChain[el]
				else :
					partname = SumSaveChain[el][0]
				indicesFull.insert(el,extindices)
				GenIndices.insert(el,partname)
				#add the external to the SummedIndices
				SummedIndices.append(extindices)
			#At this point we need to create all the factors that are in product and loop through them
			Out = []
			for ib,block in enumerate(FullTChain) : 
				if STOP :
					break
				if type(block[0]) == list : #we are dealing with a sum of terms
					#at the end all the elements of the block are summed together that is understood !!
					tempsum = []
					for iee,element in enumerate(block) : #these terms are multiplied together
						#The required number of indices is given by LinkSumSaveChain
						#copy the required indices
						prodterm = []#product of the terms inside the bracket
						#do the re ordering
						if iee != 0 : 
							iee-=1#0 has been excluded from the construction of the mapping from the very beginning
							cpGenIndices = [GenIndices[i] for i in BlockMapping[ib][iee]]+GenIndices[LinkSumSaveChain[ib]:]
							cpindicesFull = [indicesFull[i] for i in BlockMapping[ib][iee]]+ indicesFull[LinkSumSaveChain[ib]:]
						else :
							#cpGenIndices,cpindicesFull=list(GenIndices),list(indicesFull)
							cpGenIndices,cpindicesFull=copy.deepcopy(GenIndices),copy.deepcopy(indicesFull)
						for term in element :
							if term == _Th or term == _T or term == _Ta : 
								exit("The linearity for theta and or T matrices is not taken into account exiting")
							res = self.getcoeffterm(cpGenIndices,cpindicesFull,term)
							if res[0] == 0 : 
								tempsum.append(0)
								break 
							else :
								prodterm.append(res[0])
								cpGenIndices,cpindicesFull = cpGenIndices[res[-1]:],cpindicesFull[res[-1]:]
						prodterm = functools.reduce(operator.mul,prodterm) if len(prodterm) == len(element) else 0
						if prodterm != 0 :
							tempsum.append(prodterm)
					tempsum = sum(tempsum)
					if tempsum == 0 : 
						STOP = True
						break
					GenIndices,indicesFull = GenIndices[LinkSumSaveChain[ib]:],indicesFull[LinkSumSaveChain[ib]:]
					Out.append(tempsum)
				else : 
					for element in block :
						res = self.getcoeffterm(GenIndices,indicesFull,element)
						if len(res) == 3 and res[1] != [] : 
							ExtraConstrains.append(res[1])
							SummedIndices.append(res[1])
						if res[0] != 0 : 
							Out.append(res[0])
							GenIndices,indicesFull = GenIndices[res[-1]:],indicesFull[res[-1]:]
						else :
							STOP = True
							break
			if Out != [] and not(STOP): 
				#Enter the trace
				if Traces != (None,None) :
					OutOfTrace,InTrace = ([el for iel,el in enumerate(Out) 
																	if not(iel in range(Traces[0],Traces[-1]+1))],
																[el for iel,el in enumerate(Out)
																	if iel in range(Traces[0],Traces[-1]+1)])
					#get the factors out of the trace
					InTrace = functools.reduce(operator.mul,InTrace).expand()
					if type(InTrace) == Add : 
						InTS = [[arg for arg in el.args if hasattr(arg,'arg')
							 or (type(arg) == adjoint and hasattr(arg.args[0],'arg'))] for el in InTrace.args]
						#InTF = list(copy.copy(InTrace).args)
						InTF = list(copy.deepcopy(InTrace).args)
						for iel,el in enumerate(InTS):
							for ell in el : 
								InTF[iel] =InTF[iel].subs(ell,1) 
						Out = functools.reduce(operator.mul,OutOfTrace,1)*sum([
							InTF[iel]*Tr(functools.reduce(operator.mul,InTS[iel])) for iel in range(len(InTF))])
					else :
						InTS = [arg for arg in InTrace.args
								if hasattr(arg,'arg')	or (type(arg) == adjoint and hasattr(arg.args[0],'arg'))]
						#set the Symbols to one
						#InTF = list(InTrace)
						InTF = copy.deepcopy(InTrace)
						for el in InTS:
							InTF = InTF.subs(el,1)
						Out = functools.reduce(operator.mul,OutOfTrace,1)*InTF*Tr(functools.reduce(operator.mul,InTS))
				else :
					Out = functools.reduce(operator.mul,Out,1)
				listsum = [[ll for ll in el if len(str(ll[0]).split('dum')) != 2 
					and len(str(ll[1]).split('dum')) != 2] for el in SummedIndices]#got rid of the dummy indices
				#apply the extraconstrains 
				listsum = sum(listsum,[])
				constrains = [1]
				if ExtraConstrains != [] : 
					constrains = self.applyExtraConstrains2(ExtraConstrains,MatStruc)
				#multiply Out by the constrains
				Out= sum([Out * llcons for llcons in constrains])
				#create the sum out of this 
				argsum = [str(ll) for ll in listsum] 
				if argsum != [] :
					Out = eval("Sum(Out,"+','.join(argsum)+")")
				FinalTerm.append(Out)
		FinalTerm = sum(FinalTerm)
		return FinalTerm

	def Expand(self,args,dotrace=True,Layer=0,MatStruc=[]):
		"""Expand the term with valid particles
			Each term is given by ('Y',p1,p2,s1) ...
			MatStruc is used for the Mat representations. It gives the structure of the chain in this case.
			Layer counts the nested calls to Expand needed for the saved versions
			"""
		Layer += 1
		#we start by identifying al the particles that go together
		scalars = self.ScalarsN
		fermions = self.FermionsN
		groups = self.GaugeGroupsName 
		#check if there is only one argument transform the tuple
		if type(args[0]) != tuple :
			args = (args,)
		if dotrace :
			Traces = self.DeterminTraces(args)
		else : 
			Traces = (None,None)
		###################
		#end Initialization
		###################
		#Determin the structure of the chain
		SaveChain = sum([list(el[1:]) for el in args],[])
		TChain = [el[0] for el in args]
		#Check if there is a Trilinear,scalar mass, fermion mass
		TriIn = any([True if el == _h else False for el in TChain])
		FMIn = any([True if el == _mf or el == _mfa else False for el in TChain])
		ScMIn = any([True if el == _ms else False for el in TChain])
		if TriIn and self.TriToCalculate == {} :
			return 0
		if FMIn and self.FMToCalculate == {} :
			return 0
		if ScMIn and self.ScMToCalculate == {} :
			return 0
		GeneratePosition = []
		RemovedIndices = []
		PartInd = {}
		ip = 0
		#we need to identify the indices that are fixed in the chain. These appear as non global variables
		FixedIndices = [iel for iel,el in enumerate(SaveChain) if type(el) == list or not(str(el) in self.glob)]
		#remove the fix indices from the generation
		Chain = [el for iel,el in enumerate(SaveChain) if not(iel in FixedIndices)]
		#translate using the dictionary above
		TypeChain = [self.translation[el] for el in TChain]
		TypeChain = sum(TypeChain,[])
		TypeChain = [el for iel,el in enumerate(TypeChain) if not(iel in FixedIndices)]
		#do the link between the indices in the chain and the generated indices below
		for it,term in enumerate(Chain) :
			if term in PartInd:
				ind1 = PartInd[term]#index of the last identic term
				#found the corresponding value in Generate position
				ind2 = [ii[1] for ii in GeneratePosition if ii[0] == ind1]
				assert len(ind2) ==1
				GeneratePosition.append((it,ind2[0]))
				RemovedIndices.append(it)
			else :
				PartInd[term] = it#create the entry
				GeneratePosition.append((it,ip))
				ip+=1
		#create the lists for the cartesian product
		Lists = [fermions if el == 'f' else (scalars if el == 's' else groups) 
				for iel,el in enumerate(TypeChain) if not(iel in RemovedIndices)]
		#CopyLists = list(Lists)
		CopyLists = copy.deepcopy(Lists)
		#create the iterator i.e. cartesian product
		it = itr.product(*Lists)
		FinalTerm = []#for the final result
		#iteration over the combination
		for icm,combo in enumerate(it) : 
			ExtraConstrains = [] #used to do some matching for the Theta and t matrices
			STOP = False
			temp = [combo[el[1]] for el in GeneratePosition]
			Out = []
			iigen = Index(self)#create a generator of indices
			#generate the indices qui vont bien 
			indices = iigen.GenConsInd(combo,Layer)
			#duplicate the indices that are repeated
			indicesFull = [indices[el[1]] for el in GeneratePosition]#position in Generate are given with respect to the Chain list :w
			#SummedIndices = list(indices)
			SummedIndices = copy.deepcopy(indices)
			#add the ones that are fixed in the correct place
			for el in FixedIndices : 
				extindices = iigen.GenExtInd(SaveChain[el],Layer)
				if type(SaveChain[el]) == str :#group factors are not given under a list but just the name
					partname = SaveChain[el]
				else :
					partname = SaveChain[el][0]
				indicesFull.insert(el,extindices)
				temp.insert(el,partname)
				#add the external to the SummedIndices
				SummedIndices.append(extindices)
			def tempupdate(nbindtoremove,temp,indicesFull):
				#temp function which is repeated a lot
				if res != 0 :
					Out.append(res)
					temp = temp[nbindtoremove:]#remove the relevant elements
					indicesFull = indicesFull[nbindtoremove:]
					return False,temp,indicesFull
				else :
					return True,temp,indicesFull
			for itp,tp in enumerate(TChain):
				if STOP :
					break
				if tp == _Y:
					res = self.SuperY(temp[:3],indicesFull[:3],'Yukawas')
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)#update teh different dictionary and break if True is returned
				elif tp == _Ya :
					res = self.SuperY(temp[:3],indicesFull[:3],'Yukawas',Adj=True)
					res = res.adjoint() 
					STOP ,temp,indicesFull= tempupdate(3,temp,indicesFull)#update teh different dictionary and break if True is returned
				elif tp == _L : 
					res = self.SuperY(temp[:4],indicesFull[:4],'QuarticTerms')
					STOP,temp,indicesFull = tempupdate(4,temp,indicesFull)
				elif tp == _mf:
					#TODO
					res = self.SuperY(temp[:2],indicesFull[:2],'FermionMasses')
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == _mfa:
					#TODO
					res = self.SuperY(temp[:2],indicesFull[:2],'FermionMasses',Adj=True)
					res = res.adjoint()
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == _ms : 
					res = self.SuperY(temp[:2],indicesFull[:2],'ScalarMasses')
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == _h : 
					res = self.SuperY(temp[:3],indicesFull[:3],'TrilinearTerms')
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == _C or tp == _Cs : 
					res = self.Casimir(temp[:2])#we don't need the indices the casimir is indep of that
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == _S :
					res = self.DkIndex(temp[:1])#we don't need the indices the casimir is indep of that
					STOP,temp,indicesFull = tempupdate(1,temp,indicesFull)
				elif tp == _Ss :
					res = self.DkIndex(temp[:1],scalar=True)#we don't need the indices the casimir is indep of that
					STOP,temp,indicesFull = tempupdate(1,temp,indicesFull)
				elif tp == _G :
					res = self.GroupGaugeConstant(temp[:1])
					STOP,temp,indicesFull = tempupdate(1,temp,indicesFull)
				elif tp == _Cg :
					res = self.Casimir(temp[:1],Group=True)
					STOP,temp,indicesFull = tempupdate(1,temp,indicesFull)
				elif tp == _T : 
					indlabel = iigen.nextind('',Layer)
					res,Indx = self.TMatRepr(temp[:3],indicesFull[:3],indlabel)
					ExtraConstrains.append(Indx)
					if Indx != [] :
						SummedIndices.append([Indx])
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == _Ta : 
					indlabel = iigen.nextind('',Layer)
					res,Indx = self.TMatRepr(temp[:3],indicesFull[:3],indlabel,adj=True)
					ExtraConstrains.append(Indx)
					if Indx != [] :
						SummedIndices.append([Indx])
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == _Th :
					indlabel = iigen.nextind('',Layer)
					res,Indx = self.ThMatRepr(temp[:3],indicesFull[:3],indlabel)
					ExtraConstrains.append(Indx)
					if Indx != []:
						SummedIndices.append([Indx])
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == 'Y2F':
					res = self.Y2F(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'Y2Fa':
					res = self.Y2F(temp[:2],indicesFull[:2],adj=True)
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'Yab2S':
					res = self.Yab2S(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'Y2FabS':
					res = self.Y2FabS(temp[:3],indicesFull[:3])
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == 'Chain2Y':
					res = self.Chain2Yuk(temp[:4],indicesFull[:4])
					STOP,temp,indicesFull = tempupdate(4,temp,indicesFull)
				elif tp == 'Chain2Ya':
					res = self.Chain2Yuk(temp[:4],indicesFull[:4],adj=True)
					STOP,temp,indicesFull = tempupdate(4,temp,indicesFull)
				elif tp == 'Chain3Y':
					res = self.Chain3Yuk(temp[:5],indicesFull[:5])
					STOP,temp,indicesFull = tempupdate(5,temp,indicesFull)
				elif tp == 'Chain3Ya':
					res = self.Chain3Yuk(temp[:5],indicesFull[:5],adj=True)
					STOP,temp,indicesFull = tempupdate(5,temp,indicesFull)
				elif tp == 'Chain4Y':
					res = self.Chain4Yuk(temp[:6],indicesFull[:6])
					STOP,temp,indicesFull = tempupdate(6,temp,indicesFull)
				elif tp == 'Chain4Ya':
					res = self.Chain4Yuk(temp[:6],indicesFull[:6],adj=True)
					STOP,temp,indicesFull = tempupdate(6,temp,indicesFull)
				elif tp == 'Chain5Y':
					res = self.Chain5Yuk(temp[:7],indicesFull[:7])
					STOP,temp,indicesFull = tempupdate(7,temp,indicesFull)
				elif tp == 'Chain5Ya':
					res = self.Chain5Yuk(temp[:7],indicesFull[:7],adj=True)
					STOP,temp,indicesFull = tempupdate(7,temp,indicesFull)
				elif tp == 'Hbar2abS' : 
					res = self.Hbar2abS(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'H2abS' : 
					res = self.H2abS(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'Habcd':
					res = self.Habcd(temp[:4],indicesFull[:4])
					STOP,temp,indicesFull = tempupdate(4,temp,indicesFull)
				elif tp == 'Habc':
					res = self.Habc(temp[:3],indicesFull[:3])
					STOP,temp,indicesFull = tempupdate(3,temp,indicesFull)
				elif tp == 'Hab' :
					res = self.Hab(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'Theta2' : 
					res = self.Theta2(temp[:5],indicesFull[:5])
					STOP,temp,indicesFull = tempupdate(5,temp,indicesFull)
				elif tp == 'Theta4' : 
					res = self.Theta4(temp[:6],indicesFull[:6])
					STOP,temp,indicesFull = tempupdate(6,temp,indicesFull)
				elif tp == 'L2abS':
					res = self.L2abS(temp[:2],indicesFull[:2])
					STOP,temp,indicesFull = tempupdate(2,temp,indicesFull)
				elif tp == 'ffabc':
					indlabel = [iigen.nextind('',Layer) for ll in range(6)]
					res,Indx = self.fabc(temp[:1],indicesFull[:1],indlabel)
					ExtraConstrains.append(Indx)
					ExtraConstrains = sum(ExtraConstrains,[])
					if Indx != []:
						SummedIndices.append(Indx)
					STOP,temp,indicesFull = tempupdate(1,temp,indicesFull)
				else :
					exit("Expand Function : not implemented yet")
			if Out != [] and not(STOP): 
				#Enter the trace
				if Traces != (None,None) :
					OutOfTrace,InTrace = ([el for iel,el in enumerate(Out) 
																	if not(iel in range(Traces[0],Traces[-1]+1))],
																[el for iel,el in enumerate(Out)
																	if iel in range(Traces[0],Traces[-1]+1)])
					#get the factors out of the trace
					InTrace = functools.reduce(operator.mul,InTrace).expand()
					if type(InTrace) == Add : 
						InTS = [[arg for arg in el.args if hasattr(arg,'arg')
						 or (type(arg) == adjoint and hasattr(arg.args[0],'arg'))] for el in InTrace.args]
						#InTF = list(copy.copy(InTrace).args)
						InTF = list(copy.deepcopy(InTrace).args)
						for iel,el in enumerate(InTS):
							for ell in el : 
								InTF[iel] =InTF[iel].subs(ell,1) 
						Out = functools.reduce(operator.mul,OutOfTrace,1)*sum([
							InTF[iel]*Tr(functools.reduce(operator.mul,InTS[iel])) for iel in range(len(InTF))])
					else :
						InTS = [arg for arg in InTrace.args
								if hasattr(arg,'arg')	or (type(arg) == adjoint and hasattr(arg.args[0],'arg'))]
						#set the Symbols to one
						#InTF = copy.copy(InTrace)
						InTF = copy.deepcopy(InTrace)
						for el in InTS:
							InTF = InTF.subs(el,1)
						Out = functools.reduce(operator.mul,OutOfTrace,1)*InTF*Tr(functools.reduce(operator.mul,InTS))
				else :
					Out = functools.reduce(operator.mul,Out,1)
				listsum = [[ll for ll in el if len(str(ll[0]).split('dum')) != 2 
					and len(str(ll[1]).split('dum')) != 2] for el in SummedIndices]#got rid of the dummy indices
				#apply the extraconstrains 
				listsum = sum(listsum,[])
				constrains = [1]
				if ExtraConstrains != [] : 
					constrains = self.applyExtraConstrains2(ExtraConstrains,MatStruc)
				#multiply Out by the constrains
				Out= sum([Out * llcons for llcons in constrains])
				#create the sum out of this 
				argsum = [str(ll) for ll in listsum] 
				if argsum != [] :
					Out = eval("Sum(Out,"+','.join(argsum)+")")
				FinalTerm.append(Out)
		FinalTerm = sum(FinalTerm)
		return FinalTerm

	def applyExtraConstrains2(self,ExtraConstrains,MatStruc):
		"""Matches the indices """
		if type(MatStruc[0]) == str :
			MatStruc = [MatStruc]
		else : 
			pass 
		ListKronecker = []
		for ll in MatStruc : 
			Doubles = [[(ielem,iel+ielem+1) for iel,el in enumerate(ll[ielem+1:]) if el == elem ] for ielem,elem in enumerate(ll)]
			Doubles = sum((Doubles),[])
			#construct kronecker for the required indices
			GroupConstrains = [el[0].split('A')[-1] if el != [] and len(el[0].split('fabc')) ==1 else (el[0].split('fabc')[-1] if el != [] and len(el[0].split('fabc')) != 1 else 1) for el in ExtraConstrains ]
			Constrains = functools.reduce(operator.mul,[KroneckerDelta(ExtraConstrains[i][0],ExtraConstrains[j][0])
				for (i,j) in Doubles if ExtraConstrains[i] != [] and ExtraConstrains[j] != [] ],1)
			ConstrainsGroup = functools.reduce(operator.mul,[1 if GroupConstrains[i] == GroupConstrains[j]  else 0
				for (i,j) in Doubles],1)
			ListKronecker.append(Constrains*ConstrainsGroup)
		return ListKronecker

	def getcoeffterm(self,temp,indicesFull,tp) : 
		if tp == _Y:
			res = self.SuperY(temp[:3],indicesFull[:3],'Yukawas')
			return [res,3]
		elif tp == _Ya :
			res = self.SuperY(temp[:3],indicesFull[:3],'Yukawas',Adj=True)
			res = res.adjoint()
			return [res,3]
		elif tp == _L : 
			res = self.SuperY(temp[:4],indicesFull[:4],'QuarticTerms')
			return [res,4]
		elif tp == _mf:
			res = self.SuperY(temp[:2],indicesFull[:2],'FermionMasses')
			return [res,2]
		elif tp == _mfa:
			res = self.SuperY(temp[:2],indicesFull[:2],'FermionMasses',Adj=True)
			res = res.adjoint()
			return [res,2]
		elif tp == _ms : 
			res = self.SuperY(temp[:2],indicesFull[:2],'ScalarMasses')
			return [res,2]
		elif tp == _h : 
			res = self.SuperY(temp[:3],indicesFull[:3],'TrilinearTerms')
			return [res,3]
		elif tp == _C or tp == _Cs : 
			res = self.Casimir(temp[:2])#we don't need the indices the casimir is indep of that
			return [res,2]
		elif tp == _S :
			res = self.DkIndex(temp[:1])#we don't need the indices the casimir is indep of that
			return [res,1]
		elif tp == _Ss :
			res = self.DkIndex(temp[:1],scalar=True)#we don't need the indices the casimir is indep of that
			return [res,1]
		elif tp == _G :
			res = self.GroupGaugeConstant(temp[:1])
			return [res,1]
		elif tp == _Cg :
			res = self.Casimir(temp[:1],Group=True)
			return [res,1]
		elif tp == _T : 
			indlabel = iigen.nextind('',Layer)
			res,Indx = self.TMatRepr(temp[:3],indicesFull[:3],indlabel)
			return [res,Indx,3]
		elif tp == _Ta : 
			indlabel = iigen.nextind('',Layer)
			res,Indx = self.TMatRepr(temp[:3],indicesFull[:3],indlabel,adj=True)
			return [res,Indx,3]
		elif tp == _Th :
			indlabel = iigen.nextind('',Layer)
			res,Indx = self.ThMatRepr(temp[:3],indicesFull[:3],indlabel)
			return [res,Indx,3]
		elif tp == 'Y2F':
			res = self.Y2F(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'Y2Fa':
			res = self.Y2F(temp[:2],indicesFull[:2],adj=True)
			return [res,2]
		elif tp == 'Yab2S':
			res = self.Yab2S(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'Y2FabS':
			res = self.Y2FabS(temp[:3],indicesFull[:3])
			return [res,2]
		elif tp == 'Chain3Y':
			res = self.Chain3Yuk(temp[:5],indicesFull[:5])
			return [res,5]
		elif tp == 'Chain3Ya':
			res = self.Chain3Yuk(temp[:5],indicesFull[:5],adj=True)
			return [res,5]
		elif tp == 'Chain2Y':
			res = self.Chain2Yuk(temp[:4],indicesFull[:4])
			return [res,4]
		elif tp == 'Chain2Ya':
			res = self.Chain2Yuk(temp[:4],indicesFull[:4])
			return [res,4]
		elif tp == 'Hbar2abS' : 
			res = self.Hbar2abS(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'H2abS' : 
			res = self.H2abS(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'Habcd':
			res = self.Habcd(temp[:4],indicesFull[:4])
			return [res,4]
		elif tp == 'Habc':
			res = self.Habc(temp[:3],indicesFull[:3])
			return [res,3]
		elif tp == 'Hab' :
			res = self.Hab(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'Theta2' : 
			res = self.Theta2(temp[:5],indicesFull[:5])
			return [res,5]
		elif tp == 'Theta4' : 
			res = self.Theta4(temp[:6],indicesFull[:6])
			return [res,4]
		elif tp == 'L2abS':
			res = self.L2abS(temp[:2],indicesFull[:2])
			return [res,2]
		elif tp == 'ffabc':
			indlabel = [iigen.nextind('',Layer) for ll in range(6)]
			res,Indx = self.fabc(temp[:1],indicesFull[:1],indlabel)
			return [res,Indx,1]
		else :
			exit("Expand Function : not implemented yet, contact the authors")

	def DeterminTraces2(self,argsS,Structure):
		"""determin which symbols are inside a trace in the args list"""
		#get all the fermions and fermions only 
		FChain = []
		#Construct the fermion chain
		#		Structure = flatten([[ell[0] for ell in el] if type(el) == list else el for el in argsSimple])
		for ielem,elem in enumerate(argsS) :
			if elem[0] == _Y or elem[0] == _Ya :
				#get the starting value
				StartingV = tuple(elem[2]) if type(elem[2]) == list else elem[2]
			elif elem[0] == _C or elem[0] == _Cs:
				StartingV = elem[2]
			elif elem[0] == _T or elem[0] == _Ta : 
				StartingV = elem[2]
			elif elem[0] == _S or elem[0] == _Ss or elem[0] == _Cg or elem[0] == _L:
				StartingV = None
			elif elem[0] == _G or elem[0] == 'ffabc' or elem[0] == _h: #no fermion indices
				StartingV = None
			elif elem[0] == 'Y2F' or elem[0] == 'Y2Fa' or elem[0] == _mf or elem[0] == _mfa:
				StartingV = elem[1]
			elif elem[0] == 'Yab2S' or elem[0] == 'Theta2' or elem[0] == 'Theta4' or elem[0] == _Th or elem[0] == 'Agabcd':
				StartingV = None
			elif elem[0] == 'Y2FabS' or elem[0] == _ms or elem[0] == 'Hab':
				StartingV = None
			elif elem[0] == 'Hbar2abS' or elem[0] == 'Hbarabc' or elem[0] =='Habc': 
				StartingV = None
			elif elem[0] == 'H2abS' or elem[0] == 'Habcd' or elem[0] == 'L2abS' or elem[0] == 'Habc':
				StartingV = None
			elif elem[0] == 'Chain2Y' or elem[0] == 'Chain2Ya':
				StartingV = elem[3]
			elif elem[0] == 'Chain3Y' or elem[0] == 'Chain3Ya':
				StartingV = elem[4]
			elif elem[0] == 'Chain4Y' or elem[0] == 'Chain4Ya':
				StartingV = elem[5]
			elif elem[0] == 'Chain5Y' or elem[0] == 'Chain5Ya':
				StartingV = elem[6]
			else :
				exit("not implemented")
			#Get all the final values
			iel4Y = 0
			if StartingV != None :
				FinalV = []
				if len(argsS) == 1  or len(argsS[ielem:]) == 1 : 
					if argsS[ielem][0] == 'Chain4Y' or argsS[ielem][0] == 'Chain4Ya':
						iel4Y += 1
					ielem -= 1 
				for ix,x in enumerate(argsS[ielem+1:]):
					if x[0] == _Y or x[0] == _Ya : 
						FinalV.append(tuple(x[3]) if type(x[3]) == list else x[3])
					elif x[0] == _mf or x[0] == _mfa : 
						FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					elif x[0] == 'Chain4Y' or x[0] == 'Chain4Ya':
						FinalV.append(tuple(x[6]) if type(x[6]) == list else x[6])
					elif x[0] == 'Chain3Y' or x[0] == 'Chain3Ya' : 
						FinalV.append(tuple(x[4]) if type(x[5]) == list else x[5])
					elif x[0] == _mf or x[0] == _mfa or x[0] == 'Y2F' or x[0] == 'Y2Fa': 
						FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					#elif x[0] == _C :
					#	FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					#Is the initial in the Final list ? 
					if StartingV in FinalV : 
						#we found one let's get it
						InTrace= [iel for iel,el in enumerate(argsS) if iel in range(ielem+iel4Y,ix+ielem+1+1)]
						if elem[0] == _C :#we have to subtract one since C is a number multiplied out already 
							InTrace.remove(ielem)
						if elem[0] == _C : #WARNING IF THE C(F) IS OUT OF A SUMMED BLOCK IT IS NOT MULTIPLIED OUT AND THIS WOULD CRASH
							InTrace = [el - 1 for el in InTrace]
						return InTrace
		#if we get out there is no trace 
		return None,None


	def DeterminTraces(self,args):
		"""determin which symbols are inside a trace in the args list"""
		#get all the fermions and fermions only 
		FChain = []
		#Construct the fermion chain
		Structure = [el[0] for el in args]
		for ielem,elem in enumerate(args) :
			if elem[0] == _Y or elem[0] == _Ya :
				#get the starting value
				StartingV = tuple(elem[2]) if type(elem[2]) == list else elem[2]
			elif elem[0] == _C or elem[0] == _Cs:
				StartingV = elem[2]
			elif elem[0] == _T or elem[0] == _Ta : 
				StartingV = elem[2]
			elif elem[0] == _S or elem[0] == _Ss or elem[0] == _Cg or elem[0] == _L:
				StartingV = None
			elif elem[0] == _G or elem[0] == 'ffabc' or elem[0] == _h: #no fermion indices
				StartingV = None
			elif elem[0] == 'Y2F' or elem[0] == 'Y2Fa' or elem[0] == _mf or elem[0] == _mfa:
				StartingV = elem[1]
			elif elem[0] == 'Yab2S' or elem[0] == 'Theta2' or elem[0] == 'Theta4' or elem[0] == _Th or elem[0] == 'Agabcd':
				StartingV = None
			elif elem[0] == 'Y2FabS' or elem[0] == _ms or elem[0] == 'Hab':
				StartingV = None
			elif elem[0] == 'Hbar2abS' or elem[0] == 'Hbarabc' or elem[0] =='Habc': 
				StartingV = None
			elif elem[0] == 'H2abS' or elem[0] == 'Habcd' or elem[0] == 'L2abS' or elem[0] == 'Habc':
				StartingV = None
			elif elem[0] == 'Chain2Y' or elem[0] == 'Chain2Ya':
				StartingV = elem[3]
			elif elem[0] == 'Chain3Y' or elem[0] == 'Chain3Ya':
				StartingV = elem[4]
			elif elem[0] == 'Chain4Y' or elem[0] == 'Chain4Ya':
				StartingV = elem[5]
			elif elem[0] == 'Chain5Y' or elem[0] == 'Chain5Ya':
				StartingV = elem[6]
			else :
				exit("not implemented")
			#Get all the final values
			iel4Y = 0
			if StartingV != None :
				FinalV = []
				if len(args) == 1  or len(args[ielem:]) == 1 : 
					if args[ielem][0] == 'Chain4Y' or args[ielem][0] == 'Chain4Ya':
						iel4Y += 1
					ielem -= 1 
				for ix,x in enumerate(args[ielem+1:]):
					if x[0] == _Y or x[0] == _Ya : 
						FinalV.append(tuple(x[3]) if type(x[3]) == list else x[3])
					elif x[0] == _mf or x[0] == _mfa : 
						FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					elif x[0] == 'Chain4Y' or x[0] == 'Chain4Ya':
						FinalV.append(tuple(x[6]) if type(x[6]) == list else x[6])
					elif x[0] == 'Chain3Y' or x[0] == 'Chain3Ya' : 
						FinalV.append(tuple(x[4]) if type(x[5]) == list else x[5])
					elif x[0] == _mf or x[0] == _mfa or x[0] == 'Y2F' or x[0] == 'Y2Fa': 
						FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					#elif x[0] == _C :
					#	FinalV.append(tuple(x[2]) if type(x[2]) == list else x[2])
					#Is the initial in the Final list ? 
					if StartingV in FinalV : 
						#we found one let's get it
						InTrace= [iel for iel,el in enumerate(args) if iel in range(ielem+iel4Y,ix+ielem+1+1)]
						if elem[0] == _C :#we have to remove the _C index since it is transformed as a number when called it is not in the trace
							InTrace.remove(ielem)
						elif x[0] == _C : 
							inTrace.remoce(ix)
						return InTrace
		#if we get out there is no trace 
		return None,None

#####################
#End Expand functions
#####################


################################################
#Wrapping functions for Yukawas,mf,ms,Casimir...
################################################

 	def SuperY(self,parts,indices,term,Adj=False):
		"""Check if the given combination is valid and return the coeff, zero else"""
		split = False
		#check in the DB
		indices = [tuple([tuple(ll) for ll in el]) for el in indices]
		key = tuple([term] + parts + indices)
		if key in self.Databased : 
			out = self.Databased[key]
			#loggingInfo("reading superY: {}".format(term))
			if Adj == out[0][-1]:#test only the first one they should have the same
				if not(split) :
					return sum([oo[0]*oo[1] for oo in out]) 
				else :
					loggingCritical("Not implemented don't know where this is used check for compatibility with the sum",verbose=True)
					exit("")
			else : 
				if not(split) :
					Return = Integer(0)
					for oo in out : 
						if oo[0] != 0 :
							Return += oo[0].switchindices()*oo[1]
					return Return
				else :
					Return = []
					for oo in out :
						if oo[0] != 0 : 
							Return.append((oo[0].switchindices(),oo[1]))
					return Return 
		else :
			#Check that sc, f1,f2 is a valid combination
			TermIdentified = []
			if term in self.Potential : 
				for kval,val in self.Potential[term].items():
					#local copy of val
					for listfield in val :
						lcel = [str(el.args[0]) for el in listfield[0]]
						for part in parts : 
							if str(part) in lcel:
								lcel.remove(str(part)) 
							else :
								break#if one is not in the list it is useless to check the others
						if lcel == [] : 
							#the term is in this term
							TermIdentified.append((kval,listfield))
							break
			if TermIdentified == [] :
				#it is not a proper term return 0 
				#save the result
				self.Databased[key] = [(Integer(0),Integer(0),Adj)]
				return Integer(0) if not(split) else (Integer(0),Integer(0))
			else :
				#We just need to get the correct coefficient
				#When the coeff has been identified there is one more subtlety which is that the order is not conserved i.e. Pi Sig Pi Sig will be identified
				#as Sig,Sig,Pi,Pi therefore a mapping is needed
				ToDatabased = []
				ToOut = []
				for tterm in TermIdentified : 
					lcparts = [str(elem) for elem in parts]
					NewInd = [elem for elem in indices]
					lctterm = [str(elem.args[0]) for elem in tterm[1][0]]
					while True :
						NotAllowed = []
						for ix,x in enumerate(lcparts) : #go through the not ordered list and find the first place where it should be 
							if ix in NotAllowed :
								pass
							else :
								for new,nn in enumerate(lctterm) :
									if new in NotAllowed :
										pass
									else :
										if x == nn and not(new in NotAllowed) :
											if ix == new : 
												NotAllowed.append(ix)
												break
											else :
												if lcparts[ix] != lcparts[new] :
													lcparts[ix] = lcparts[new]
													lcparts[new] = x
													buff = NewInd[ix]
													NewInd[ix] = NewInd[new]
													NewInd[new] = buff
													NotAllowed.append(ix)
													NotAllowed.append(new)
													break
												else :
													break
						if lctterm == lcparts:
							break
						else :
							pass
					#collect the different symbols for each groups
					FinalOrdering = [[el[il] for el in NewInd] for il in range(len(self.NonUGaugeGroups))] 
					#get rid of the dummy indices
					FinalIndices = [[ll[0] for ll in el if len(str(ll[0]).split('dum')) != 2] for el in FinalOrdering]
					#Get rid of the singlets that are charge under no group and that therefore appears as empty lists
					FinalIndices = [el for el in FinalIndices if el != []]
					#construct the factor
					if type(tterm[1][-1]) != int :
						Factor = tterm[1][-1].subs(sum([[(tterm[1][1][iid][iie][0],ee) for iie,ee in enumerate(ii)]for iid,ii in enumerate(FinalIndices)],[]))
					else :
						Factor = tterm[1][-1]#Singlet case the factor is not a list but an int
					if term == 'Yukawas' or term =='FermionMasses':
						#This is just to know the position of the femrions in the parts list
						i0 = 1 if term == 'Yukawas' else 0
						if Adj : 
							ToDatabased.append((self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0+1],parts[i0]),Factor,Adj))
						else :
							ToDatabased.append((self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0],parts[i0+1]),Factor,Adj))
					else :
						ToDatabased.append((Symbol(tterm[0],commutative=True),Factor,Adj))
					if term == 'Yukawas' or term == 'FermionMasses': 
						i0 = 1 if term == 'Yukawas' else 0
						if Adj : 
							if split :
								ToOut.append((self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0+1],parts[i0]),Factor)) 
							else :
								ToOut.append((self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0+1],parts[i0])*Factor))
						else :
							if split :
								ToOut.append(self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0],parts[i0+1]),Factor)
							else :
								ToOut.append(self.Classes[tterm[0]](Symbol(tterm[0],commutative=False),parts[i0],parts[i0+1])*Factor)
					else :
						if split :
							ToOut.append((Symbol(tterm[0],commutative=True),Factor))
						else :
							ToOut.append(Symbol(tterm[0],commutative=True)*Factor)
				self.Databased[key] = ToDatabased
				out = sum(ToOut)
			return out


	def Casimir(self,label,Group=False):
		"""wrapper for the casimir operator."""
		#the particle and teh group are always known at this point the Casimir can be calculated
		key = tuple(['C2(F)']+label) if not(Group) else tuple(['C2(G)']+label)
		if key in self.InvariantResults:
			return self.InvariantResults[key]
		else :
			if not(Group) :
				res = self.GetGroupFromName[label[0]].C2(self.Particles[str(label[1])].Qnb[label[0]])
			else :
				res = self.GetGroupFromName[label[0]].C2('G')
			self.InvariantResults[key] = res
			return res

	def DkIndex(self,label,scalar=False):
		"""Calculate the Dynkin Index for fermions or scalars under the given gaugegroup"""
		key = tuple(['S(F)']+label) if not(scalar) else tuple(['S(S)']+label)
		if key in self.InvariantResults:
			return self.InvariantResults[key]
		else :
			if not(scalar) :
				res = DynkinIndex(self,[label[0],self.GetGroupFromName[label[0]]],'Fermions')
			else :
				res = DynkinIndex(self,[label[0],self.GetGroupFromName[label[0]]],'Scalars')
			self.InvariantResults[key] = res
			return res

	def GroupGaugeConstant(self,gp):
		"""return the g value that is correct"""
		return self.GetGroupFromName[gp[0]].g**2


	def fabc(self,label,indices,indlabel):
		gr = self.GetGroupFromName[label[0]]
		#create a symbolic Function
		if not(gr.U):
			fabc = gr.fabc
			indlabel = [['{}fabc{}'.format(el,gr.name),0,gr.d-1] for el in indlabel]
			res = Fabc(fabc[0],fabc[1],indlabel[0][0],indlabel[1][0],indlabel[2][0],indlabel[3][0],indlabel[4][0],indlabel[5][0])
			return res,indlabel
		else :
			return 0,[]


	def ThMatRepr(self,label,indices,indlabel,adj=False):
		"""wrap the theta matrices and calculates them """
		key = tuple(['ThMat']+label+flatten(indices)+[indlabel]) if not(adj) else tuple(['ThaMat']+label+flatten(indices)+[indlabel])
		if key in self.InvariantResults:
			return self.InvariantResults[key]
		else :
			#the two scalars have to come from the same higgs particle 
			if type(label[1]) == str :
				label[1] = Symbol(label[1])
			if type(label[2]) == str :
				label[2] = Symbol(label[2])
			Higgs1 = set(self.getHiggsFromScalar[label[1]])
			Higgs2 = set(self.getHiggsFromScalar[label[2]])
			Inter = Higgs1.intersection(Higgs2) != set([])
			grp = self.GetGroupFromName[label[0]]
			if Inter : 
				if not(grp.U) :
					for idgr,(nm,gr,uu) in enumerate(self.NonUGaugeGroups):
						if nm == grp.name : 
							break	
					#charge under the gauge group
					f1 = [self.ScalarsIndex[label[1]]]+[self.Particles[str(label[1])].Qnb[label[0]]]+ [ind for idi,ind in enumerate(indices[1]) if idi == idgr]
					f2 = [self.ScalarsIndex[label[2]]]+[self.Particles[str(label[2])].Qnb[label[0]]] + [ind for idi,ind in enumerate(indices[2]) if idi == idgr]
					if len(str(f1[-1][0]).split('dum')) == 2 : #particles are not charged under the group
						self.InvariantResults[key] = Integer(0),[]
						return Integer(0),[]
					else :
						#Get the other groups indices which have to be conserved
						f1kr = [ind for idi,ind in enumerate(indices[1]) if idi != idgr]
						f2kr = [ind for idi,ind in enumerate(indices[2]) if idi != idgr]
						#create a symbolic function for the matrix
						#generate the index 
						Indx = [('{}A{}'.format(indlabel,label[0])),0,grp.d-1]
						#create the product of kronecker to conserve the other qnb
						Ind = [[f1kr[i][0],f2kr[i][0]] for i in range(len(f1kr)) if len(str(f1kr[i][0]).split('dum')) != 2]
						Kr = functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Ind],1)
						#Now for the scalars we use the HB basis i.e. if the two scalars are the same they are multiplied by mat-mat*
						#if they are different they are multiplied by +-(mat+mat*) 
						res = Ts(Indx[0],grp.MatRepr(f1,f2,HB=True),f1[-1][0],f2[-1][0])*Kr  #take the indices label for f1 and f2
						res = res if not(adj) else res.conjugate()
						self.InvariantResults[key] = res,Indx
						return res,Indx
				else :#U1 group
					f1,f2 = label[1],label[-1]
					f1ind = [el[0] for el in indices[1] if len(str(el[0]).split('dum')) != 2]
					f2ind = [el[0] for el in indices[2] if len(str(el[0]).split('dum')) != 2]
					#Get the Hypercharge
					f1h = self.Particles[str(label[1])].Qnb[label[0]]#f2h is the same condition enforced at the begining
					f1,f2 = [self.ScalarsIndex[label[1]],f1h] + f1ind,[self.ScalarsIndex[label[2]],f1h] + f2ind
					res = grp.MatRepr(f1,f2,HB=True)#no need for symbolic function here everything is known and nothing summed easy !
					res = res if not(adj) else res.conjugate()
					self.InvariantResults[key] = res,[]
					return res,[]
			else :
				self.InvariantResults[key] = Integer(0),[]
				return Integer(0),[]


	def TMatRepr(self,label,indices,indlabel,adj=False):
		"""Wrap the T matrices and calculates them"""
		#TODO To make it faster one can get teh ind out of the function but be carefull doing so check that the generated indices are diffrent
		key = tuple(['TMat']+label+flatten(indices)+[indlabel]) if not(adj) else tuple(['TaMat']+label+flatten(indices)+[indlabel])
		if key in self.InvariantResults:
			return self.InvariantResults[key]
		else :
			if label[1] == label[-1] :#They have to come from the same particle
				grp = self.GetGroupFromName[label[0]]
				if not(grp.U) :
					for idgr,(nm,gr,uu) in enumerate(self.NonUGaugeGroups):
						if nm == grp.name : 
							break	
					#check that the two particles are in from the same particle and that they are charge under the gauge group
					f1 = [self.Particles[str(label[1])].Qnb[label[0]]]+ [ind for idi,ind in enumerate(indices[1]) if idi == idgr]
					f2 = [self.Particles[str(label[2])].Qnb[label[0]]] + [ind for idi,ind in enumerate(indices[2]) if idi == idgr]
					#Get the other groups indices which hqve to be conserved
					f1kr = [ind for idi,ind in enumerate(indices[1]) if idi != idgr]
					f2kr = [ind for idi,ind in enumerate(indices[2]) if idi != idgr]
					if len(str(f1[-1][0]).split('dum')) == 2 : #particles are not charged under the group
						self.InvariantResults[key] = Integer(0),[]
						return Integer(0),[]
					else :
						#create a symbolic function for the matrix
						#generate the index 
						Indx = [('{}A{}'.format(indlabel,label[0])),0,grp.d-1]
						#create the product of kronecker to conserve the other qnb
						Ind = [[f1kr[i][0],f2kr[i][0]] for i in range(len(f1kr)) if len(str(f1kr[i][0]).split('dum')) != 2]
						Kr = functools.reduce(operator.mul,[KroneckerDelta(i,j) for (i,j) in Ind],1)
						res = Ts(Indx[0],grp.MatRepr(f1[0],f2[0]),f1[-1][0],f2[-1][0])*Kr  #take the indices label for f1 and f2
						res = res if not(adj) else res.conjugate()
						self.InvariantResults[key] = res,Indx
						return res,Indx
				else :#U1 group
					f1,f2 = label[1],label[-1]
					f1ind = [el[0] for el in indices[1] if len(str(el[0]).split('dum')) != 2]
					f2ind = [el[0] for el in indices[2] if len(str(el[0]).split('dum')) != 2]
					#Get the Hypercharge
					f1h = self.Particles[str(label[1])].Qnb[label[0]]#f2h is the same condition enforced at the begining
					f1,f2 = [f1h] + f1ind,[f1h] + f2ind
					res = grp.MatRepr(f1,f2)#no need for symbolic function here everything is known and nothing summed easy !
					res = res if not(adj) else res.conjugate()
					self.InvariantResults[key] = res,[]
					return res,[]
			else :
				self.InvariantResults[key] = Integer(0),[]
				return Integer(0),[]


		
##########################################
#DEFINITION of the INVARIANTS
##########################################
	
	def getparts(self,parts,indices):
		lenparts = len(parts)
		Pindices = [[el[0] for el in indices[i]] for i in range(lenparts)]
		Parts = [[el] + Pindices[iel] for iel,el in enumerate(parts)]
		return Parts

	def Y2F(self,parts,indices,adj=False):
		"""Calculates the invariant Y2F"""
		key = tuple(['Y2F']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading Y2F")
		else :
			f1,f2 = self.getparts(parts,indices)
			if not(adj):
				res = self.Expand(((_Ya,s1,f1,p1),(_Y,s1,p1,f2)),Layer=1,dotrace=False)
			else :
				res = self.Expand(((_Y,s1,f1,p1),(_Ya,s1,p1,f2)),Layer=1,dotrace=False)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Yab2S(self,parts,indices):
		"""calculates the invariant Yab2S Eq. 19"""
		key = tuple(['Yab2S']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading Yab2S")
		else :
			sc1,sc2 = self.getparts(parts,indices)
			res = Rational(1,2)*(
					self.Expand(((_Ya,sc1,p1,p2),(_Y,sc2,p2,p1)),Layer=1)
				+ self.Expand(((_Ya,sc2,p1,p2),(_Y,sc1,p2,p1)),Layer=1) 
				)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Habc(self,parts,indices):
		"""Calculates the invariant Habc Eq 66"""
		key = tuple(['Habc']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading Habcd")
		else :
			sc1,sc2,sc3 = self.getparts(parts,indices)
			res =	[Rational(1,2)*(self.Expand(((_mf,p1,p2),(_Ya,a,p2,p3),(_Y,b,p3,p4),(_Ya,c,p4,p1)),Layer=1)
							+ self.Expand(((_Y,a,p1,p2),(_mfa,p2,p3),(_Y,b,p3,p4),(_Ya,c,p4,p1)),Layer=1)
							)
							for (a,b,c) in list(permutations([sc1,sc2,sc3],3))
							]
			res = sum(res)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Habcd(self,parts,indices):
		"""Calculates the invariant Habcd Eq 40"""
		key = tuple(['Habcd']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading Habcd")
		else :
			sc1,sc2,sc3,sc4 = self.getparts(parts,indices)
			res = [Rational(1,4)*(
					self.Expand((('Chain4Y',a,b,c,d,p1,p1)))
					)
					for (a,b,c,d) in list(permutations([sc1,sc2,sc3,sc4],4))
					]
			res = sum(res)
			if res.subs(Tr(0),0) != 0 : 
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Hab(self,parts,indices):
		"""Calculates the invariant Eq 87"""
		key = tuple(['Hab']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading Hab")
		else :
			sc1,sc2= self.getparts(parts,indices)
			#res = (self.Expand(((_Y,sc1,p1,p2),(_Ya,sc2,p2,p3),(_mf,p3,p4),(_mfa,p4,p1)),Layer=1)
			#		+ self.Expand(((_Y,sc2,p1,p2),(_Ya,sc1,p2,p3),(_mf,p3,p4),(_mfa,p4,p1)),Layer=1)
			#		+ self.Expand(((_Ya,sc1,p1,p2),(_Y,sc2,p2,p3),(_mfa,p3,p4),(_mf,p4,p1)),Layer=1)
			#		+ self.Expand(((_Ya,sc2,p1,p2),(_Y,sc1,p2,p3),(_mfa,p3,p4),(_mf,p4,p1)),Layer=1)
			#		+ self.Expand(((_Y,sc1,p1,p2),(_mfa,p2,p3),(_Y,sc2,p3,p4),(_mfa,p4,p1)),Layer=1)
			#		+ self.Expand(((_mf,p1,p2),(_Ya,sc1,p2,p3),(_mf,p3,p4),(_Ya,sc2,p4,p1)),Layer=1)
			#		)
			res = (self.Expand2((([[(_Y,sc1,p1,p2),(_Ya,sc2,p2,p3)],[(_Y,sc2,p1,p2),(_Ya,sc1,p2,p3)]]),(_mf,p3,p4),(_mfa,p4,p1)),Layer=1)#Line 1
					+ self.Expand2((([[(_Ya,sc1,p1,p2),(_Y,sc2,p2,p3)],[(_Ya,sc2,p1,p2),(_Y,sc1,p2,p3)]]),(_mfa,p3,p4),(_mf,p4,p1)),Layer=1)#Line 2
					+ self.Expand(((_Y,sc1,p1,p2),(_mfa,p2,p3),(_Y,sc2,p3,p4),(_mfa,p4,p1)),Layer=1)
					+ self.Expand(((_mf,p1,p2),(_Ya,sc1,p2,p3),(_mf,p3,p4),(_Ya,sc2,p4,p1)),Layer=1)
					)
			if res != 0 : 
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	
	def L2abS(self,parts,indices):
		"""Calculates Eq 24"""
		key = tuple(['L2abS']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading L2abs")
		else :
			a,b = self.getparts(parts,indices)
			res = Rational(1,6)*self.Expand(((_L,a,s1,s2,s3),(_L,b,s1,s2,s3)),Layer=1)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Y2FabS(self,parts,indices):
		"""Calculates the invariant Y2FabS Eq. 27"""
		key = tuple(['Y2FabS']+parts+flatten(indices))
		if key in self.InvariantResults:
			loggingDebug("reading Y2FabS")
			res = self.InvariantResults[key]
		else :
			grp,sc1,sc2 = self.getparts(parts,indices)
			res = Rational(1,2)*(
					self.Expand(((_C,grp,p1),(_Ya,sc1,p1,p2),(_Y,sc2,p2,p1)),Layer=1)
				+ self.Expand(((_C,grp,p1),(_Ya,sc2,p1,p2),(_Y,sc1,p2,p1)),Layer=1) 
				)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res
	
	def Theta4(self,parts,indices) : 
		key = tuple(['Theta4']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else :
			g1,g2,sc1,sc2,sc3,sc4 = self.getparts(parts,indices)
			res = ( self.Expand(((_Th,g1,sc1,s1),(_Th,g1,sc3,s2),(_Th,g2,s1,sc2),(_Th,g2,s2,sc4)),Layer=1,MatStruc=['A','A','B','B'])
				+ self.Expand(((_Th,g1,sc1,s2),(_Th,g1,s1,sc4),(_Th,g2,s1,sc3),(_Th,g2,sc2,s2)),Layer=1,MatStruc=['A','A','B','B']))
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res


	def Theta2(self,parts,indices):#not used 
		key = tuple(['Theta2']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else : 
			g1,g2,sc1,sc2 = self.getparts(parts,indices)
			res = self.Expand(((_Th,g1,sc1,s1),(_Th,g2,s1,sc2)),Layer=1,MatStruc=['A','B'])
			if res != 0 : 
				res = res.doit()
				self.InvariantResults[key] = res
		return res



	def Chain2Yuk(self,parts,indices,adj=False):
		key = tuple(['Chain2Y']+parts+flatten(indices)) if not(adj) else tuple(['Chain2Ya']+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else :
			sc1,sc2,f1,f2 = self.getparts(parts,indices)
			if not(adj) : 
				res = self.Expand(((_Y,sc1,f1,p2),(_Y,sc2,p2,f2)),Layer = 1, dotrace=False)
			else :
				res = self.Expand(((_Ya,sc1,f1,p1),(_Y,sc2,p1,f2)),Layer = 1, dotrace=False)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	
	def Chain3Yuk(self,parts,indices,adj=False):
		"""Calculates a Chain of Three yukawas"""
		key = tuple(['Chain3Y']+parts+flatten(indices)) if not(adj) else tuple(['Chain3Ya']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
			loggingDebug("reading chain3Y")
		else :
			sc1,sc2,sc3,f1,f2 = self.getparts(parts,indices)
			if not(adj) : 
				res = self.Expand2(((_Y,sc1,f1,p1),(_Ya,sc2,p1,p2),(_Y,sc3,p2,f2)),Layer = 1, dotrace=False)
			else :
				res = self.Expand2(((_Ya,sc1,f1,p1),(_Y,sc2,p1,p2),(_Ya,sc3,p2,f2)),Layer = 1, dotrace=False)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Chain4Yuk(self,parts,indices,adj=False):
		"""Calculates a Chain of four Yukawas"""
		key = tuple(['Chain4Y']+parts+flatten(indices)) if not(adj) else tuple(['Chain4Ya']+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else :
			sc1,sc2,sc3,sc4,f1,f2 = self.getparts(parts,indices)
			if not(adj) : 
				res = self.Expand((('Chain3Y',sc1,sc2,sc3,f1,p1),(_Ya,sc4,p1,f2)),Layer=2,dotrace=False)
			else :
				res = self.Expand((('Chain3Ya',sc1,sc2,sc3,f1,p1),(_Y,sc4,p1,f2)),Layer=2,dotrace=False)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res

	def Chain5Yuk(self,parts,indices,adj=False):
		"""Calculates a Chain of five Yukawas"""
		key = tuple(['Chain5Y']+parts+flatten(indices)) if not(adj) else tuple(['Chain5Ya']+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else :
			sc1,sc2,sc3,sc4,sc5 = parts[:5]
			f1,f2 = parts[5:]
			sc1ind = [el[0] for el in indices[0]]
			sc2ind = [el[0] for el in indices[1]]
			sc3ind = [el[0] for el in indices[2]]
			sc4ind = [el[0] for el in indices[3]]
			sc5ind = [el[0] for el in indices[4]]
			f1ind = [el[0] for el in indices[5]]
			f2ind = [el[0] for el in indices[6]]
			sc1,sc2,sc3,sc4,sc5 = [sc1]+sc1ind,[sc2]+sc2ind,[sc3]+sc3ind,[sc4]+sc4ind,[sc5]+sc5ind
			f1,f2 = [f1]+f1ind,[f2]+f2ind
			if not(adj) : 
				res = self.Expand((('Chain4Y',sc1,sc2,sc3,sc4,f1,p1),(_Y,sc5,p1,f2)),Layer=3,dotrace=False)
			else :
				res = self.Expand((('Chain4Ya',sc1,sc2,sc3,sc4,f1,p1),(_Ya,sc5,p1,f2)),Layer=3,dotrace=False)
			if res != 0 :
				res = res.doit()
			self.InvariantResults[key] = res
		return res


	def Hbar2abS(self,parts,indices):
		"""Calculates Eq 26 """
		key = tuple(['Hbar2abS']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else : 
			sc1,sc2 = parts
			sc1ind = [el[0] for el in indices[0]]
			sc2ind = [el[0] for el in indices[1]]
			sc1,sc2 = [sc1]+sc1ind,[sc2]+sc2ind
			res = Rational(1,2)*(
					self.Expand((('Chain3Y',sc1,s2,sc2,p1,p2),(_Ya,s2,p2,p1)),Layer=1)
				+	self.Expand((('Chain3Ya',sc1,s2,sc2,p1,p2),(_Y,s2,p2,p1)),Layer=1)
				)
			if res != 0 :
				res = res.doit()	
			self.InvariantResults[key] = res
		return res

	def H2abS(self,parts,indices):
		"""Calculates Eq 25 """
		key = tuple(['H2abS']+parts+flatten(indices))
		if key in self.InvariantResults:
			res = self.InvariantResults[key]
		else : 
			sc1,sc2 = parts
			sc1ind = [el[0] for el in indices[0]]
			sc2ind = [el[0] for el in indices[1]]
			sc1,sc2 = [sc1]+sc1ind,[sc2]+sc2ind
			res = Rational(1,2)*(
					self.Expand2((('Chain3Y',sc1,sc2,s2,p1,p2),(_Ya,s2,p2,p1)),Layer=1)
				+ self.Expand2((('Chain3Ya',sc1,sc2,s2,p1,p2),(_Y,s2,p2,p1)),Layer=1)
				)
			if res != 0 :
				res = res.doit()	
			self.InvariantResults[key] = res
		return res

####################################################
#END Wrapping functions for Yukawas,mf,ms,Casimir...
####################################################


	def declareSymbol(self,string) : 
		"""declare globally the symbol contained in string"""
		#add the symbol to the global dictionary
		if type(string) == list : 
			for st in string : 
				#modifies the pyrate module globals() dictionary
				self.glob[st] = Symbol(st)
				#this modify the module globals() dictionary 
				globals()[st] = Symbol(st)
		else :
			#modifies the pyrate module globals() dictionary
			self.glob[string] = Symbol(string)
			#this modify the module globals() dictionary 
			globals()[string] = Symbol(string)

	def declarations(self):
			"""does all the global declarations"""
			#particle place holders
			self.declareSymbol(flatten([('p{}'.format(i),'s{}'.format(i),'gg{}'.format(i)) for i in range(15)]))
			#Yukawas,Casimir,Dynkin,Group,mf,ms,Lambda,trilinear,thethas,tAs
			self.declareSymbol(['_Y','_Ya','_C','_Cs','_S','_Ss','_Cg','_G','_mf','_mfa','_ms','_L','_h','_Tha','_Ta','_Th','_T']) #the last symbls is for the C2(G)
			#declare all the yukawas and quartics appearing in the model file
			if 'Yukawas' in self.Potential :
				self.declareSymbol(self.Potential['Yukawas'].keys())

	def CheckSpecifiedParticlesScalars(self,listParticles):
		"""Check if the given list of particles contains only defined particles"""
		#Check if they are in the HiggsFields
		if all([el in self.Potential['HiggsFields'] for el in listParticles]) : 
			return True
		else : 
			return all([el in self.Scalars for el in listParticles])

###################
#Indices Generation
###################


class Index(object):
	"""deals with the indices"""

	def __init__(self,model):
		self.GeneratedIndices = []
		self.mod = model
		self.intind = 0
		self.extind = 0
		self.grintind = 0
		self.grextind = 0
		self.dummy = 0
	
	def __repr__(self) : 
		string = "Index instance :\n\t\t\tinternal {}\n\t\t\texternal {}\n\t\t\tgroup".format(self.intind,self.extind,self.grintind)
		return string

	def GenConsInd(self,combo,Layer):
		#generate a set of consistent indices
		out = []
		for ip,pp in enumerate(combo):
			lencombo = len(combo)
			Gout = []
			if pp in self.mod.GaugeGroupsName:
				Gout.append([self.nextind(pp,Layer,group=True),pp,pp])
			else :
				for group in self.mod.NonUGaugeGroups:
					#get teh particle
					ip1 = self.mod.Particles[str(pp)].rangeindices[group[0]]
					if ip1 != 0 :#singlet otherwise
						Gout.append([self.nextind(pp,Layer),1,ip1])
					else :
						Gout.append([self.nextind(pp,Layer,dummy=True)])
			out.append(Gout)
		return out

	def GenExtInd(self,part,Layer):
		"""this take a particle and given isospin indices and create an ext indices that is not sum but only a single value"""
		Gout = []
		#deal with the group factors
		if type(part) == str and part in self.mod.GaugeGroupsName:
			Gout.append([self.nextind(part,Layer,ext=True,group=True),part,part])
		elif type(part) != str and part[0] in self.mod.GaugeGroupsName:
			Gout.append([self.nextind(part[0],Layer,ext=True,group=True),part[1],part[1]])
		else :
			for ig,group in enumerate(self.mod.NonUGaugeGroups):
				if part[ig+1] != 0 or type(part[ig+1]) == Symbol :#the particle is charged
					if len(str(part[ig+1]).split('dum'))== 2 :#check if the ext ind is map to a dummy one in this case keep dumm
						Gout.append([self.nextind(part[0],Layer,dummy=True),part[ig+1],part[ig+1]])
					else :
						Gout.append([self.nextind(part[0],Layer,ext=True),part[ig+1],part[ig+1]])
				else :#dummy
					Gout.append([self.nextind(part[0],Layer,ext=True,dummy=True),part[ig+1],part[ig+1]])
		return Gout

	def nextind(self,pp,Layer,ext=False,dummy=False,group=False) :
		#exterior indices are not sum in the first place
		if not(ext):
			if not(group) and not(dummy):
				self.intind += 1 
				tempstring = 'int{}{}{}'.format(self.intind,Layer,str(pp))
				self.mod.declareSymbol(tempstring)
				return tempstring
			elif dummy : #for particles not charged under some groups
				self.dummy += 1
				tempstring = 'dum{}{}'.format(self.dummy,Layer)
				self.mod.declareSymbol(tempstring)
				return tempstring
			else :
				self.grintind += 1
				tempstring = 'intgg{}{}{}'.format(self.grintind,Layer,pp)
				self.mod.declareSymbol(tempstring)
				return tempstring
		else : 
			if not(dummy) and not(group): 
				self.extind +=1
				tempstring ='ext{}{}{}'.format(self.extind,Layer,str(pp))
				self.mod.declareSymbol(tempstring)
				return tempstring
			elif not(dummy) and group : 
				self.grextind +=1
				tempstring = 'exgg{}{}{}'.format(self.grextind,Layer,pp)
				self.mod.declareSymbol(tempstring)
				return tempstring
			else : 
				self.dummy +=1
				tempstring = 'dum{}{}'.format(self.dummy,Layer)
				self.mod.declareSymbol(tempstring)
				return tempstring

	def reset(self) :
		self.intind = 0
		self.extind = 0
		self.dummy = 0 
		self.grextind = 0 
		self.grintind = 0 
				


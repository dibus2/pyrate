#!/usr/bin/env python
from __future__ import division
import numpy as np
from numpy import sqrt,log,exp
import scipy.integrate as scint
import scipy.interpolate as sci
import scipy.optimize as sco

"""date: 9th May 2013	I define a RGE class in order to store the RGEs and then solve them."""

class PhysicalParameters(object) : 
	"""This class contains all the different physical parameters for the calculations but also some parameters for the run such as the higgs mass.	All the constants are given at the Z  boson mass."""
	def __init__(self):
		"""Set the physical constants to their values."""
		self.Pi = np.pi
		self.alpha = 1./128.91
		self.alphaS = 0.1184
		self.sinthetasq = 0.2316 #Sin squared
		self.Mz = 91.1876 #in GeV
		self.Mw = 80.399 # in GeV
		self.Gf = 1.16637e-5 #in GeV ^-2
		self.CF = 4./3. #Casimir of the SU3 color
		self.Nc = 3.
		self.Qtop = 2./3.#top quark chargs
		self.vev = 246.22 # From (sqrt(2)Gf)^-1/2
		self.Mt = 173.3 # the pole mass of the top
		self.U = [2.4e-3,1.27,165.3308]#mu,mc,mt (MSbar)mt=(Mt/(1+4alphas/3Pi))
		self.D = [4.75e-3,104e-3,4.19]#md,ms,mb
		self.L = [0.511e-3,105.66e-3,1.777]
		self.yU = [sqrt(2)/self.vev * el for el in self.U]
		self.yD = [sqrt(2)/self.vev * el for el in self.D]
		self.yL = [sqrt(2)/self.vev * el for el in self.L] 
		self.g10 = sqrt(4*self.Pi * self.alpha)/sqrt(1-self.sinthetasq)#*sqrt(5./3.)#SU5 norm#TODO
		self.g20 = sqrt(4*self.Pi * self.alpha)/sqrt(self.sinthetasq)
		self.g30 =  sqrt(4*self.Pi * self.alphaS)

	def __repr__(self):
		"""Represent the class instances."""
		return "Instance of the RGE class. Physical parameters predefined calculated at the Z boson mass : \n\t Pi {}\n\t alpha {}\n\t alphaS {}\n\t sinthetasq {} \n\t ZbosonMass {}\n\t Gf {}\n\t Higgs vev {}\n\t Top pole mass {}\n\t Up matrix {}\n\t Down matrix {}\n\t Lepton matrix {}\n.".format(self.Pi,self.alpha,self.alphaS,self.sinthetasq,self.Mz,self.Gf,self.vev,self.Mt,self.U,self.D,self.L)

class RGE(PhysicalParameters) :
	"""inHerits from PhysicalParameters all the parameters."""
	"""This class contains the RGEs and the method to solve them."""

	def __init__(self,beta,nbeqs,labels=[],Y0=[]):
		"""Initialize the RGE.Takes Y0 and the beta function"""
		PhysicalParameters.__init__(self)
		#Solutions of the RGE
		self.Y0 = Y0
		self.nbeqs = nbeqs
		self.beta = beta
		self.assumptions = {'two-loop': True, 'diag': True}
		self.Sol = {}
		self.Solint = {}#interpolated results
		if labels == [] :
			self.labels = range(self.nbeqs)
		else :
			assert len(labels) == self.nbeqs
			self.labels = labels

	def betafunction(self,t,y):
		"""wrapping of the beta function in order to pass the assunption argument, impossible to get it work whith the set_param method"""
		if self.assumptions != {} :
			return self.beta(t,y,Assumptions = self.assumptions)
		else :
			return self.beta(t,y)
		
	def solve_rges(self,t0,tmax,step,assumptions={},ode='dopri5'):
		"""Solve ODE using the general class for ODE solver of scipy."""
		r = scint.ode(self.betafunction).set_integrator('{}'.format(ode),atol='1e-8')
		if assumptions != {}:
			self.assumptions = assumptions
		#Do transformation of the t values into log(10) 
		#t=log(Q/Mz) tp = log10(Q)
		t0 = np.log(10)*t0 - np.log(self.Mz)
		tmax = np.log(10)*tmax - np.log(self.Mz)
		r.set_initial_value(self.Y0,t0)
		#we declare a dic to store the results
		for ii in self.labels : 
			self.Sol[ii] = []
		t = []
		while r.successful() and r.t < tmax :
			r.integrate(r.t+step)
			for idi,ii in enumerate(r.y) : 
				self.Sol[self.labels[idi]].append(ii)
			t.append((r.t+np.log(self.Mz))/np.log(10))
		#at the end we store the results	
		self.Sol['t'] = np.array(t)
		for ii in self.labels : 
			self.Sol[ii] = np.array(self.Sol[ii])
			self.Solint[ii] = sci.interp1d(self.Sol['t'],self.Sol[ii])

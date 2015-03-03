# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <rawcell>

# 
# This is a script that calls mathematica and builds up the database from Susyno results.

# <codecell>

import string
import pickle
import re
import copy
import subprocess
import itertools,functools
import pudb
# <codecell>

#open the file as a string do the replacement un put it back before loading it
f = open('CGCsnew.m','r')
dic = f.readlines()
f.close()

# <codecell>

dic = ''.join(dic).replace('^','**')
while 'Sqrt[' in dic:
    dic =dic.replace('Sqrt[{}]'.format(re.search('Sqrt\[([^\]]*)\]',dic).group(1)),'sqrt({})'.format(re.search('Sqrt\[([^\]]*)\]',dic).group(1)))
while re.search('sqrt\([0-9]*/[0-9]*\)',dic) != None :
    src = re.search('sqrt\(([0-9]*)/([0-9]*)\)',dic)  
    dic = dic.replace('sqrt({}/{})'.format(src.group(1),src.group(2)),'sqrt(Rational({},{}))'.format(src.group(1),src.group(2)))#
while re.search('[0-9]{2}/[0-9]{2}',dic) != None :
   src = re.search('([0-9]{2})/([0-9]{2})',dic)
   dic = dic.replace('{}/{}'.format(src.group(1),src.group(2)),'(Rational({},{}))'.format(src.group(1),src.group(2)))
while re.search('[0-9]{1}/[0-9]{2}',dic) != None :
   src = re.search('([0-9]{1})/([0-9]{2})',dic)
   dic = dic.replace('{}/{}'.format(src.group(1),src.group(2)),'(Rational({},{}))'.format(src.group(1),src.group(2)))
while re.search('[0-9]{2}/[0-9]',dic) != None :
   src = re.search('([0-9]{2})/([0-9])',dic)
   dic = dic.replace('{}/{}'.format(src.group(1),src.group(2)),'(Rational({},{}))'.format(src.group(1),src.group(2)))
while re.search('[0-9]/[0-9]',dic) != None :
   src = re.search('([0-9])/([0-9])',dic)
   dic = dic.replace('{}/{}'.format(src.group(1),src.group(2)),'(Rational({},{}))'.format(src.group(1),src.group(2)))
# <codecell>

f = open('CGCsnew.py','w')
f.write('from sympy import sqrt,Rational,I\n')
f.write(dic)
f.close()
# <codecell>

from sympy import sqrt,Rational,I,Matrix
from CGCsnew import db

# <codecell>

#Now we can play around with the databased and polish the resul
#What has to be done is the following : 1.1 
tp = {}
#Correct SU2
#for key,val in db['SU2']['Bilinear'].items():
#    kk = ((key[0],),(key[1],)),((key[0],True),(key[1],)), ((key[0],),(key[1],True)), ((key[0],True),(key[1],True))
#    tp[kk[0]],tp[kk[-1]] = val,val
#    tp[kk[1]],tp[kk[2]] = [tuple([i]*2+[1]) for i in range(1,len(val)+1)],[tuple([i]*2+[1]) for i in range(1,len(val)+1)]
#db['SU2']['Bilinear'] = copy.deepcopy(tp)
#Construct the Quartic terms :
#for group,val in db.items() :
#	if group != 'date' and group != 'Matrices' and group != 'DynkinToDim' and group != 'Dynkin' and group != 'Casimir':
#		Quartickeys = list(itertools.product(*[[list(el) for el in val['Bilinear'].keys()],[list(el) for el in val['Bilinear'].keys()]]))
#		Quartickeys = [tuple(el[0] + el[1]) for el in Quartickeys]
#		Quarticelem = list(itertools.product(*[val['Bilinear'].values(),val['Bilinear'].values()]))
#		for iel,el in enumerate(Quarticelem) :
#			TP = []
#			prod = list(itertools.product(*[el[0],el[1]]))
#			for iell,ell in enumerate(prod) :
#				tp = [0]*(len(ell[0])+len(ell[1]) - 1)
#				tp[:len(ell[0])-1] = ell[0][:-1]
#				tp[len(ell[0])-1: len(ell[1])+len(ell[0])-2] = ell[1][:-1]
#				tp[-1] = ell[0][-1] * ell[1][-1]
#				TP.append(tuple(tp))
#			Quarticelem[iel] = copy.deepcopy(TP)
#		db[group]['Quartic'] = {}
#		for ikk,kk in enumerate(Quartickeys):
#			db[group]['Quartic'][kk] = Quarticelem[ikk] 
tp = {}
tpbar = {}
for key,val in db['DynkinToDim']['SU2'].items():
    tp[(key,)] = val
    tp[(key,True)] = val
    tpbar[val],tpbar[-val] = (key,),(key,True)
db['SU2']['DimToDynkin'] = copy.deepcopy(tpbar)
db['SU2']['DimToDynkin'][1] = (0,)#singlet
db['SU2']['DynkinToDim'] = copy.deepcopy(tp)
db['SU2']['DynkinToDim'][(0,)] = 0
#All other groups: We just need to construct the DimToDynkin dict
for key,val in db['DynkinToDim'].items():
    tp = {}
    if key != 'SU2':
        for kk,vv in val.items():
            tp[vv] = kk
            if kk[::-1] != kk:
                tp[-vv]=kk[::-1]#Assume symmetry (0,0,0,0,1)-->(1,0,0,0)
        db[key]['DimToDynkin'] = copy.deepcopy(tp)
        db[key]['DimToDynkin'][1] = tuple([0]*(int(key.split('SU')[-1])-1))
        db[key]['DynkinToDim'] = copy.deepcopy(val)
        db[key]['DynkinToDim'][tuple([0]*(int(key.split('SU')[-1])-1))] = 0
pudb.set_trace()
for key,val in db['Matrices'].items() :
    tp = {}
    tpA = {}
    if key == 'SU2' :
        db['SU2']['Fond'] = (1,)
        db['SU2']['Adj'] = (2,)
        for kk,vv in val.items():
            tp[(kk,)] = {'shape':(len(vv[0][0]),len(vv[0])),'mat': vv} 
            tpA[(kk,True)] = {'shape':(len(vv[0][0]),len(vv[0])),'mat': [tuple([tuple(-Matrix(el).conjugate().row(i)) for i in range(len(Matrix(el).col(0)))]) for el in vv]} 
    else :
        fond = [0]*(int(key.split('U')[-1])-1)
        fond[0] = 1
        adj = copy.deepcopy(fond)
        adj[-1] = 1
        fond = tuple(fond)
        db[key]['Fond'] = fond
        db[key]['Adj'] = tuple(adj)
        
        for kk,vv in val.items():
            tp[kk] ={'shape': (len(vv[0][0]),len(vv[0])), 'mat': vv}
            #tpA[kk] = {'shape':(len(vv[0][0]),len(vv[0])),'mat': [tuple([tuple(-Matrix(el).conjugate().row(i)) for i in range(len(Matrix(el).col(0)))]) for el in vv]} 
    for kk,val in tp.items() :
      tp[kk]['mat'] = tuple(tp[kk]['mat'])
    for kk,val in tpA.items() :
      tpA[kk]['mat'] = tuple(tpA[kk]['mat'])
    db[key]['Matrices'] = copy.deepcopy(tp) 
    db[key]['Matrices'].update(tpA) 


for key,val in db['Casimir'].items():
	db[key]['Casimir'] = copy.deepcopy(val)
for key,val in db['Dynkin'].items():
	db[key]['Dynkin'] = copy.deepcopy(val)


# <codecell>
del(db['DynkinToDim'])
del(db['Matrices'])
del(db['Casimir'])
del(db['Dynkin'])


# <codecell>


#Let's now calculate the HB basis as well as the AMatrices for the irreps in the models
from sympy import Matrix,zeros
for group,val in db.items():
    if group != 'date' :
        val['HBmat'] = {}
        val['HBAmat'] = {}
        for irreps,vval in val['Matrices'].items() :
            Tp,TpA = [],[]
            val['HBmat'][irreps] = {}
            val['HBAmat'][irreps] = {}
            tp, tpA = {},{}
            shp = vval['shape']
            tp[(0,0)] = {'mat': [], 'shape' :shp}
            tpA[(0,0)] = {'mat': [], 'shape' :shp}
            tp[(0,1)] = {'mat': [], 'shape' :shp}
            tpA[(0,1)] = {'mat': [], 'shape' :shp}
            tp[(1,0)] = {'mat': [], 'shape' :shp}
            tpA[(1,0)] = {'mat': [], 'shape' :shp}
            tp[(1,1)] = {'mat': [], 'shape' :shp}
            tpA[(1,1)] = {'mat': [], 'shape' :shp}
            for mat in vval['mat'] :
                #Construct the matrix
                tp[(0,0)]['mat'].append(Rational(1,2)*(Matrix(mat)-Matrix(mat).conjugate()))
                tp[(0,1)]['mat'].append(Rational(1,2)*I*(Matrix(mat)+Matrix(mat).conjugate()))
                tp[(1,0)]['mat'].append(-Rational(1,2)*I*(Matrix(mat)+Matrix(mat).conjugate()))
                tp[(1,1)]['mat'].append(Rational(1,2)*(Matrix(mat)-Matrix(mat).conjugate()))
                for akey,aval in tp.items() :
                    #tpA[(akey[0],True)]['mat'].append(tuple(-aval['mat'][-1].conjugate()))
                  tp[akey]['mat'][-1] = tuple(aval['mat'][-1])
            tp[(0,0)]['mat'],tp[(0,1)]['mat'],tp[(1,0)]['mat'],tp[(1,1)]['mat'] = tuple(tp[(0,0)]['mat']),tuple(tp[(0,1)]['mat']),tuple(tp[(1,0)]['mat']),tuple(tp[(1,1)]['mat'])
            #reorganize the matrices we want in [(0,0)] all the matrices
            val['HBmat'][irreps] = copy.deepcopy(tp)
# <codecell>
#To sparse matrices
NewMatrix = {}
for group,val in db.items():
    if group != 'date' :
			NewMatrix[group] = {}
			for irreps,vval in val['Matrices'].items():
				AllMat = []
				for matrix in vval['mat']:
					tempMat = []
					for line in range(vval['shape'][0]):
						for col in range(vval['shape'][1]):
							if matrix[line][col] != 0 :
								tempMat.append((line,col,matrix[line][col]))
					AllMat.append(tuple(tempMat))
				NewMatrix[group][irreps] = copy.deepcopy(tuple(AllMat))

import numpy as np
NewHBmatrix = {}
NewAHBmatrix = {}
for group,val in db.items():
    if group != 'date' :
			NewHBmatrix[group] = {}
			print("group: {}".format(group))
			for irreps,vval in val['HBmat'].items():
				print("\tirrep translation: {}".format(irreps))
				NewHBmatrix[group][irreps] = {}
				for struc,vstruc in vval.items():
					print("\t\tstruc: {}".format(struc))
					AllMat = []
					for matrix in vstruc['mat']:
						tempMat = []
						for line in range(vstruc['shape'][0]):
							for col in range(vstruc['shape'][1]):
								tpmatrix = np.array(matrix).reshape((vstruc['shape'][0],vstruc['shape'][1]))
								if tpmatrix[line][col] != 0 :
									tempMat.append((line,col,tpmatrix[line][col]))
						AllMat.append(tuple(tempMat))
						print("\t\t\t matrix +1...")
					NewHBmatrix[group][irreps][struc] = copy.deepcopy(tuple(AllMat))

for group,val in db.items():
    if group != 'date':
			NewAHBmatrix[group] = {}
			print("group: {}".format(group))
			for irreps,vval in val['HBAmat'].items():
				print("\tirrep translation: {}".format(irreps))
				NewAHBmatrix[group][irreps] = {}
				for struc,vstruc in vval.items():
					print("\t\tstruc: {}".format(struc))
					AllMat = []
					for matrix in vstruc['mat']:
						tempMat = []
						for line in range(vstruc['shape'][0]):
							for col in range(vstruc['shape'][1]):
								tpmatrix = np.array(matrix).reshape((vstruc['shape'][0],vstruc['shape'][1]))
								if tpmatrix[line][col] != 0 :
									tempMat.append((line,col,tpmatrix[line][col]))
						AllMat.append(tuple(tempMat))
						print("\t\t\t matrix +1...")
					NewAHBmatrix[group][irreps][struc] = copy.deepcopy(tuple(AllMat))


#Calculate the structure functions for each group, old method GetStructurefunctions of Model object
from sympy import *
StructureFunctions = {}
def sumperso(list) : 
		out = list[0]
		for el in list[1:]:
			out += el
		return out

check = True
structures = {}
for group, val in db.items():
	if group != 'date':
		print "Group: {}".format(group)
		N = int(group.split('SU')[-1])
		d = N**2-1
		W = MatrixSymbol('W',d,1)
		V = MatrixSymbol('V',d,1)
		Vec1 = sumperso([W[i,0]*Matrix(val['Matrices'][val['Fond']]['mat'][i]) for i in range(d)])
		Vec2 = sumperso([V[i,0]*Matrix(val['Matrices'][val['Fond']]['mat'][i]) for i in range(d)])
		print"\t calculating ResTrace..."
		ResTrace =[sum([
			-2*(Matrix(val['Matrices'][val['Fond']]['mat'][i])[j,k]*KroneckerDelta(l,m) - Matrix(val['Matrices'][val['Fond']]['mat'][i])[l,m]*KroneckerDelta(k,j))*Vec1[k,l]*Vec2[m,j]
			for j in range(N)
			for k in range(N)
			for l in range(N)
			for m in range(N)
			if (Matrix(val['Matrices'][val['Fond']]['mat'][i])[j,k]*KroneckerDelta(l,m) - Matrix(val['Matrices'][val['Fond']]['mat'][i])[l,m]*KroneckerDelta(k,j))*Vec1[k,l]*Vec2[m,j] != 0])
			for i in range(d)]
		print"\t done."
		Structures = [[[I*(el.diff(W[i,0])).diff(V[j,0]) for i in range(d)] for j in range(d)]for el in ResTrace]
		Structures = np.array(Structures).reshape(d,d,d)
		print"\t Structure finished, check..."
		structures[group] = copy.deepcopy(Structures)
		if check :
			Check = [Matrix(val['Matrices'][val['Fond']]['mat'][i])*Matrix(val['Matrices'][val['Fond']]['mat'][j]) - Matrix(val['Matrices'][val['Fond']]['mat'][j])*Matrix(val['Matrices'][val['Fond']]['mat'][i]) + sumperso([I*Structures[i,j,k]*Matrix(val['Matrices'][val['Fond']]['mat'][k])
					for k in range(d)])
					for i in range(d)
					for j in range(d)
					]
			Check = all([el == zeros((N,N)) for el in Check])
			if Check:
				print"Check ok"
			else : 
				print "ERROR while determining the structure constants"
		else :
			pass

#I comment the sparse matrices for the moment thats not what is important
#for group,val in db.items():
#	if group != 'date':
#		del(db[group]['Matrices'])
#		del(db[group]['HBmat'])
#		del(db[group]['HBAmat'])
#		db[group]['Matrices'] = copy.deepcopy(NewMatrix[group])
#		db[group]['HBmat'] = copy.deepcopy(NewHBmatrix[group])
#		db[group]['HBAmat'] = copy.deepcopy(NewAHBmatrix[group])
		db[group]['Struc'] = copy.deepcopy(structures[group])




subprocess.call(["rm","CGCsnew.py"])

 #<codecell>

f = open('CGCsnew.pickle','w')
pickle.dump(db,f)
f.close()

# <codecell>

g = open('CGCsnew.pickle','r')
data = pickle.load(g)
g.close()

# <codecell>



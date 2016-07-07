"""
Combines different group database into a single one, e.g. su2 su3 su4 ...
"""

import os
import sys 
import pickle
import pudb
import gzip
sys.path.append(os.path.abspath('../../src/pylie/'))
from PyLie import *
from sympy import symbols, Symbol, Rational, sqrt, IndexedBase, Matrix, \
    Wild, Symbol, Function, symbols, pi, Rational, zeros, I, sqrt, eye, MatrixSymbol, \
    KroneckerDelta, flatten, pprint, IndexedBase, Idx, Integer, Add, Mul, Indexed, Sum, conjugate, adjoint, \
    __version__, Mod

from sympy import MutableMatrix as tMM
from sympy.physics.secondquant import evaluate_deltas
#######
# Shadowing of the piecewise mechanism of sympy that makes the code way too slow
import sympy.concrete.expr_with_limits
import sympy.concrete.delta

argv = sys.argv
if len(argv) <=2 :
    exit("At least two arguments needed")
filenames = argv[1:]
pickles = [pickle.load(file(ff,'r')) for ff in filenames[1:]]
with gzip.open(filenames[0], 'rb') as fdb:
    output = pickle.load(fdb)
fdb.close()
pudb.set_trace()
for ip,pp in enumerate(pickles): 
    print("adding {} to the database {}".format(filenames[ip+1], filenames[0]))
    for key,val in pp.items():
        if key in output:
            if type(val) == dict:
                for kkey,vval in val.items():
                    if kkey in output[key]:
                        if type(vval) == dict:
                            for kkkey, vvval in vval.items():
                                output[key][kkey][kkkey] = vvval
                        else:
                            output[key][kkey] = vval
                    else:
                        output[key][kkey] = vval
            else:
                output[key] = val
        else:
            output[key] = val

with gzip.open('CGCslast.pklz','wb') as fdb:            
    pickle.dump(output,fdb)
print("done")

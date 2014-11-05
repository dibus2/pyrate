# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

#This is a patch for version 0.7.3 of sympy 
#Warning this is modifying version 0.7.3 of sympy see sympy0.7.3.patch for details on the modifications 

# <codecell>

import stat
import os
import sys
# <codecell>

if 'dry' in sys.argv :
    dry = True
else :
    dry = False
if 'R' in sys.argv :
    reverse = True
else :
    reverse = False

# <codecell>

try : 
    import sympy as sp
except ImportError as err: 
    exit("error your sympy is not installed")

# <codecell>

#Get the installation path of the sympy package
path = sp.__path__
version = sp.__version__

#check version
if version != '0.7.3' :
	exit("your Sympy version is not 0.7.3, you don't need to apply the patch !")
print("checking sympy version...ok")

#check permissions for the two files that we want to modify

st1 = path[0]+'/matrices/expressions/matexpr.py'
st2 = path[0]+'/concrete/delta.py'
try :
    f = open(st1,'a')
except :
    exit("permission denied, try to run sudo python patch.sympy-0.7.3.py")
f.close()
try :
    f = open(st2,'a')
except :
    exit("permission denied, try to run sudo python patch.sympy-0.7.3.py")



# <codecell>

if dry and not(reverse):
    command_line = "patch --dry-run -p1 -d {} <sympy-0.7.3.patch".format(path[0])
elif dry and reverse :
    command_line = "patch --dry-run -R -p1 -d {} <sympy-0.7.3.patch".format(path[0])
elif not(dry) and reverse :
    command_line = "patch -R -p1 -d {} <sympy-0.7.3.patch".format(path[0])
else :
    command_line = "patch -p1 -d {} <sympy-0.7.3.patch".format(path[0])
try :
	a = os.popen(command_line)
	print(''.join(a.readlines()))
	print("Sympy-0.7.3 patched for PyR@TE !")
except :
	print("Error while patching Sympy-0.7.3")

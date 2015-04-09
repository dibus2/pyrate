"""
Combines different group database into a single one, e.g. su2 su3 su4 ...
"""



import sys 
import pickle
import pudb

argv = sys.argv
if len(argv) <=2 :
    exit("At least two arguments needed")
filenames = argv[1:]
output = {}
pickles = [pickle.load(file(ff,'r')) for ff in filenames]
for ip,pp in enumerate(pickles) : 
    print("adding {} to the database".format(filenames[ip]))
    output.update(pp)

pickle.dump(output,file('CGCslast.pickle','w'))
print("done")

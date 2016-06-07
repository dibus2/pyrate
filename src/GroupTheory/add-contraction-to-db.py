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
pickles = [pickle.load(file(ff,'r')) for ff in filenames]
output = pickles[0]
for ip,pp in enumerate(pickles[1:]): 
    print("adding {} to the database {}".format(filenames[ip+1], filenames[0]))
    for key,val in pp.items():
        pudb.set_trace()
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

            
pickle.dump(output,file('CGCslast.pickle','w'))
print("done")

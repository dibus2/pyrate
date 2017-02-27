"""
This is the RunPyLie class that allows to calculate all the group theoretic information without using the interactive mode.
"""
import sys
import os

wd = os.getcwd()
sys.path.append(wd + '/src/Core')
from IPyrate import *


class RunPyLie:
    def __init__(self, filename, output=''):
        self.filename = filename
        self._idb = Idbquerry()  # creates an instance of the Idbquerry class
        self._idb.noprint=True  # switch off the printing
        self._queries = []  # to store the list of queries
        self._lines = []  # to store the lines from the file
        self._comments = []
        self._results = []  # to store the results from the different querries
        self._output = output  # the name of the output filename
        ####
        # Execute the run
        ####
        self._setUp()  # set up the calculation

    def _loadQueries(self):
        # load the queries from the file
        if os.path.isfile(self.filename):
            with open(self.filename, 'r') as f:
                self._lines = f.readlines()
                f.close()
        else:
            loggingCritical("Error in `RunPyLie` class, file `{}` does not exist.".format(self.filename), verbose=True)
            exit()

    def _filterQueries(self):
        # remove possible lines that are comments and then makes sure that the key words are known to PyLie before passing them even though I think it would raise errors
        self._queries = [(el.split(' ')[0],el.split(' ')[1:]) for el in self._lines if el[0] != '#' and el != '\n']
        self._comments = [el for el in self._lines if el[0] == '#']
        self._lines = [el for el in self._lines if el[0] != '#' and el != '\n']


    def _setUp(self):
        # sets up the calculation and in particular reads the queries from file
        self._loadQueries()  # load the queries from the file
        self._filterQueries()  # filters the queries

    def calcAll(self):
        # goes through the queries and sends them to PyLie
        for qr in self._queries:
            res = eval("self._idb.do_{}('{}')".format(qr[0], self._idb.toline(qr[1])[:-1]))
            if not res:
                qr[1].insert(0, qr[0])
                loggingCritical("Error, calculating term `{}`.".format(' '.join(qr[1])))
            self._results.append(res)
        # need to filter the results and write them into a nice form
        for irr, rr in enumerate(self._results):
            if isinstance(rr, basestring):
                rr = [el[1:] for el in rr.split('\n') if el != '' and not 'invariant' in el]
                self._results[irr] = rr
        # export the results to txt
        if self._output == '':
            if '.' in self.filename:
                ext = self.filename.split('.')[-1]
                self._output = '.'.join(self.filename.split('.')[:-1]) + '.output.txt'
            else:
                self._output = self.filename + '.ouput.txt'
        with open(self._output, 'w') as fout:
            for iel, el in enumerate(self._lines):
                fout.write('# '+el)
                if isinstance(self._results[iel], list):
                    fout.write('invs{} = {{'.format(iel+1))
                    fout.write(','.join(self._results[iel]))
                    fout.write('}\n')
                else:
                    if self._results[iel]:
                        pudb.set_trace()
                        fout.write(self._queries[iel][0]+'{} = '.format(iel) + str(self._results[iel])+'\n')
                    else:
                        fout.write('invalid entry\n')
            fout.close()
        self._idb._update_db()  # update the database since we are here.





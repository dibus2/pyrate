#!/usr/bin/env python 
"""loading of a store expression"""
try:
    import string
    import sys

    sys.path.append('../GroupTheory')
except:
    exit("error while loading modules in ToMathematica.py")

from RGEsmathModule import *


def findclosingbracket(lstring, pos):
    # initialize counter of open brackets
    counter = 0
    while True:
        if lstring[pos] == ')' and counter == 0:
            break
        elif lstring[pos] == '(':
            counter += 1
        elif lstring[pos] == ')':
            counter -= 1
        if pos == len(lstring) - 1 and counter != 0:
            raise Exception('unbalance round brackets')
        pos += 1
    return pos


def extractSymbol(func, inStr):
    # find the position of the opening bracket
    pos0 = inStr.find('{}('.format(func)) + len(func)
    lstring = list(inStr)
    # replace by square bracket
    lstring[pos0] = '['
    # increase pos to skip the squared bracket that we have just replaced in the loop below
    pos = pos0 + 1
    # and now find the matching round bracket
    pos = findclosingbracket(lstring, pos)
    # replace the last round bracket
    lstring[pos] = ']'
    # At this stage we have to check if we are not in the case with no indices
    # Get the position of the next coma
    inStr = ''.join(lstring)
    return inStr


def extractTuple(func, inStr):
    # find the position of the opening bracket
    pos = inStr.find('{}[('.format(func)) + len(func) + 1
    # make a list of chars, so that we can make positional assignments
    lstring = list(inStr)
    # replace by square bracket
    lstring[pos] = ''
    # increase pos to skip the squared bracket that we have just replaced in the loop below
    pos += 1
    # and now find the matching round bracket
    pos = findclosingbracket(lstring, pos)
    # replace the last round bracket,we close the square bracket
    lstring[pos] = ']'
    if lstring[pos + 1] == ']':  # then there re no indices and we stop
        lstring[pos + 1] = ''
    else:
        lstring[pos + 1] = '['  # replace the coma by a squared bracket
    return string.join(lstring, '')


def roundToSquareForFunction(func, inStr):
    # by T. Jezo :)
    # find the position of the opening bracket
    pos = inStr.find('{}('.format(func)) + len(func)
    # make a list of chars, so that we can make positional assignments
    lstring = list(inStr)
    # replace by square bracket
    lstring[pos] = '['
    # increase pos to skip the squared bracket that we have just replaced in the loop below
    pos += 1
    # and now find the matching round bracket
    pos = findclosingbracket(lstring, pos)
    # replace the last round bracket
    lstring[pos] = ']'
    return string.join(lstring, '')


def ToMathematicaNotation(vall, model, FlagSquare=True):
    vall = str(vall)
    if FlagSquare:
        while 'conjugate(' in vall:
            vall = roundToSquareForFunction('conjugate', vall)
        while 'MatM(' in vall:
            vall = roundToSquareForFunction('MatM', vall)
        while 'SP(' in vall:
            vall = roundToSquareForFunction('SP', vall)
        while 'Dagger(' in vall:
            vall = roundToSquareForFunction('Dagger', vall)
        while 'trace(' in vall:
            vall = roundToSquareForFunction('trace', vall)
        while 'transpose(' in vall:
            vall = roundToSquareForFunction('transpose', vall)
        while 'sqrt(' in vall:
            vall = roundToSquareForFunction('sqrt', vall)
        while 'Sqrt(' in vall:
            vall = roundToSquareForFunction('Sqrt', vall)
    while 'MatM[(' in vall:
        vall = extractTuple('MatM', vall)
    for yuksymb in model.ListYukawa:
        while '{}('.format(yuksymb) in vall:
            vall = extractSymbol(yuksymb, vall)
    for fmsymb in model.ListFM:
        while '{}('.format(fmsymb) in vall:
            vall = extractSymbol(fmsymb, vall)
    vall = vall.replace('conjugate', 'conj').replace('MatM', 'MatMul').replace('SP', 'ScalarProd').replace('Dagger',
                                                                                                           'Adj').replace(
        'transpose', 'Tp').replace('[ ]', '')
    # Simplify the Symbols to get rid of the LateX structure
    symbs = model.ListYukawa + model.ListFM + model.ListLbd + model.ListScM + model.ListTri
    ll = []
    for el in symbs:
        if len(reg.split('{(.*)}', el)) == 3 or len(reg.split('\\\(.*)', el)) == 3:
            ll.append(''.join(reg.split('\\\(.*)', ''.join(reg.split('{(.*)}', el)))))
        else:
            ll.append(el)
    ll = [el.replace('_', '').replace('\\', '').replace('\\\\', '').replace('^','').replace(' ', '') for el in ll]
    Replace = [(el, ll[iel]) for iel, el in enumerate(symbs)]
    for rep in Replace:
        vall = vall.replace(*rep)
    return vall


def TranslateToMathematica(xpr, fname, tag, model):
    Final = []
    # returned for the AllMathematicaRGEs
    Return, Dimension = [], []
    for ll, vall in xpr.items():
        key = ll
        vall = ToMathematicaNotation(vall, model)
        if len(reg.split('{(.*)}', str(ll))) == 3 or len(reg.split('\\\(.*)', str(ll))) == 3:
            lltp = ''.join(reg.split('\\\(.*)', ''.join(reg.split('{(.*)}', str(ll))).replace('_', '')))
            ll = lltp.lower()
            ll = ll[0].upper() + ll[1:]
        else:
            ll = str(ll)
            # notsure what this is doing...
            lltp = ll
            ll = ll[0].upper() + ll[1:]
        ll = ll.replace('_', '').replace('\\', '').replace('\\\\', '').replace(' ', '')
        Final.append([ll + '=', str(vall).replace('**', '^').replace('_', '').replace('sqrt', 'Sqrt')])
        if ll in model.GetGroupFromName:
            Return.append("{{{},{}}}".format(
                str(model.GetGroupFromName[ll].g).replace('_', '').replace('\\', '').replace('**', '^'),
                str(vall).replace('**', '^').replace('_', '')))
            Dimension.append("{{{},{{}}}}".format(
                str(model.GetGroupFromName[ll].g).replace('_', '').replace('\\', '').replace('**', '^')))
        elif key in model.YukToCalculate.keys() + model.FMToCalculate.keys():
            dims = model.DimYuk[key]
            if len(dims) == 1:  # vector
                dims = [dims[0], 1]
            elif dims == []:
                dims = [1, 1]
            if dims != [1, 1]:
                elem = model.YukToCalculate[key][0] if key in model.YukToCalculate else model.FMToCalculate[key][0]
                fermions = [str(el.args[0]) + 'f' for el in elem if
                            str(el.args[0]) in model.Fermions and model.Fermions[str(el.args[0])].Gen != 1]
                if len(fermions) >= 2 and fermions[0] == fermions[1]:
                    fermions[1] = fermions[1].split('f')[0] + '1f'
                Return.append("{{{}[{}],{}}}".format(lltp.replace('_', '').replace('\\', ''), ','.join(fermions),
                                                     str(vall).replace('**', '^').replace('_', '').replace('\\', '')))
            else:
                Return.append("{{{},{}}}".format(lltp.replace('_', '').replace('\\', ''),
                                                 str(vall).replace('**', '^').replace('_', '').replace('\\', '')))
            if dims != [1, 1]:
                Dimension.append("{{{},{{{},{}}}}}".format(lltp.replace('_', '').replace('\\', ''), dims[0], dims[1]))
            else:
                Dimension.append("{{{},{{}}}}".format(lltp.replace('_', '').replace('\\', '')))
        else:  # lambdas and so on
            Return.append("{{{},{}}}".format(lltp.replace('_', '').replace('\\', '').replace('^',''),
                                             str(vall).replace('**', '^').replace('_', '').replace('\\', '')))
            Dimension.append("{{{},{{}}}}".format(lltp.replace('_', '').replace('\\', '').replace('^','')))
    f = open('{}'.format(fname), 'w')
    date = time.localtime()
    date = "{}-{}-{}\t {}:{}:{}".format(date[0], date[1], date[2], date[3], date[4], date[5])
    f.write(
        "######################################################################################\n#This is the generated by pyR@TE on {}, results for the {} terms\n######################################################################################\n\n\n\n".format(
            date, tag))
    for el in Final:
        f.write(el[0] + '\t' + el[1] + '\n')
    f.close()
    return Return, Dimension

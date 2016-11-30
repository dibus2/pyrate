#!/usr/bin/env python 
from RGEsDefinition import *
from ModelsClass import Model, loggingInfo, loggingDebug, loggingCritical
from RGEsmathModule import DynkinIndex, DynkinCasimir, SimplifyTraces, DeterminOrdering, Tr


###########################
# Fermion Anomalous dimension
###########################


def set_globalsFermionAnomalous(model):
    """coppy the globals from RGEsModule into the local gloabls()"""
    for key,val in model.glob.items():
        globals()[key] = val

def CompileFermionAnomalous(Translated,mqt,comb,model,Weyl,GutNorm):
    LH = model.FermionAnomalousToCalculate[mqt][1]
    FinalBeta = (sum([el[0] for el in Translated[mqt]]).doit()/LH).expand()
    if Weyl :
        FinalBeta = FinalBeta.subs(kappa,Rational(1,2))
    else :
        FinalBeta = FinalBeta.subs(kappa,1)
    if GutNorm:  # Apply any GUT normalization
        FinalBeta = FinalBeta.subs(model.GutNorm)
    FinalBeta = FinalBeta.subs(tuple([(el,0) for el in ListAllSymbols['FermionAnomalous']]))
    FinalBeta = DeterminOrdering(model,FinalBeta)
    return FinalBeta

def CFermionAnoI(powe,comb,model):
    """Calculates Eq 28"""
    assert powe[1] == 0
    f1,f2 = comb
    res = Rational(1,2)*model.Expand(((_Y,s1,f1,p2),(_Ya,s1,p2,f2)))
    powe[0] = powe[0].subs(FermionAnoI,res)
    return powe[0]

def CFermionAnoIg2(powe,comb,model):
    """Calculates Eq 28"""
    assert powe[1] == 2
    f1,f2 = comb
    conserv = [KroneckerDelta(a,b) for a,b in zip(f1,f2)]
    if type(conserv[0]) == KroneckerDelta  :
            conserv[0] = 0 if conserv[0].args[0] != conserv[0].args[1] else 1
    conserv = functools.reduce(operator.mul,conserv)
    if conserv != 0 :
        res = model.Expand(((_G,gg1),(_C,gg1,f1)))*conserv
    else :
        res = 0
    powe[0] = powe[0].subs(FermionAnoIg2,res)
    return powe[0]

#######
#2-loop
#######



def CFermionAnoII(powe,comb,model):
    """Calculates Eq. 29 line 1"""
    assert powe[1] == 0
    f1,f2 = comb
    res = (-Rational(1,8)*model.Expand2(((_Y,s1,f1,p1),(_Ya,s2,p1,p2),(_Y,s2,p2,p3),(_Ya,s1,p3,f2)),dotrace=False)
            -kappa*Rational(3,2)*model.Expand2(((_Y,s1,f1,p1),(_Ya,s2,p1,f2),('Yab2S',s1,s2)),dotrace=False)
            )
    if res != 0 : 
        res = res.doit()
    powe[0] = powe[0].subs(FermionAnoII,res)
    return powe[0]


def CFermionAnoIIg2(powe,comb,model):
    """Calculates Eq. 29 line 2 and 3"""
    assert powe[1] == 2
    f1,f2 = comb
    res = (Rational(9,2)*model.Expand2(((_G,gg1),(_Cs,gg1,s1),(_Y,s1,f1,p1),(_Ya,s1,p1,f2))) #line 2
            -Rational(7,4)*model.Expand2(((_G,gg1),(_C,gg1,f1),(_Y,s1,f1,p1),(_Ya,s1,p1,f2))) #line 2
            -Rational(1,4)*model.Expand2(((_G,gg1),(_Y,s1,f1,p1),(_C,gg1,p1),(_Ya,s1,p1,f2))) #line 3
            )
    if res !=0 :
        res = res.doit()
    powe[0] = powe[0].subs(FermionAnoIIg2,res)
    return powe[0]



def CFermionAnoIIg4(powe,comb,model):
    """Calculates Eq. 29 line 4 5 """
    assert powe[1] == 4 
    f1,f2 = comb
    conserv = [KroneckerDelta(a,b) for a,b in zip(f1,f2)]
    if type(conserv[0]) == KroneckerDelta  :
                    conserv[0] = 0 if conserv[0].args[0] != conserv[0].args[1] else 1
    conserv = functools.reduce(operator.mul,conserv)
    if conserv != 0 :
        res = ((Rational(25,4)+Rational(1,4)*xi**2 + 2*xi)*model.Expand2(((_G,gg1),(_G,gg1),(_C,gg1,f1),(_Cg,gg1)))
               - 2*kappa*model.Expand2(((_G,gg1),(_G,gg1),(_C,gg1,f1),(_S,gg1)))
               -Rational(1,4)*model.Expand2(((_G,gg1),(_G,gg1),(_C,gg1,f1),(_Ss,gg1)))
               -Rational(3,2)*model.Expand2(((_G,gg1),(_G,gg2),(_C,gg1,f1),(_C,gg2,f1)))
               )
    else :
        res = 0
    
    if res != 0 : 
        res = res.doit()
    powe[0] = powe[0].subs(FermionAnoIIg4,res)
    return powe[0]


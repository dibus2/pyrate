# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

def filtertraceMatM(xpr,dic):
    """
    check if there is a zero matrix in dic and return zero if it is the case.
    """
    if any([val==0 for key,val in dic.items()]) :
        return Integer(0)
    
def SimplifytraceMatM(xpr) :
    """
    remove all the traces and MatM, SP that contain a zero matrix in.
    """
    #terms with two or four objects
    p,q,r,s,t,u,v,w = map(Wild,['p','q','r','s','t','u','v','w'])
    xpr = xpr.replace(lambda xpr : xpr.match(trace(p,q)) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(trace(p,q))))
    xpr = xpr.replace(lambda xpr : xpr.match(trace(p,q,r,s)) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(trace(p,q,r,s))))
    xpr = xpr.replace(lambda xpr : xpr.match(MatM((p,q))) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(MatM((p,q)))))
    xpr = xpr.replace(lambda xpr : xpr.match(MatM((p,q,r,s))) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(MatM((p,q,r,s)))))
    xpr = xpr.replace(lambda xpr : xpr.match(MatM((p,q),r,s)) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(MatM((p,q),r,s))))
    xpr = xpr.replace(lambda xpr : xpr.match(MatM((p,q,u,v),r,s)) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(MatM((p,q,u,v),r,s))))
    xpr = xpr.replace(lambda xpr : xpr.match(SP(p,q)) != None, lambda xpr : filtertraceMatM(xpr,xpr.match(SP(p,q))))    
    return xpr

def settozero(xpr,symbs):
    """
    put symbols to zeros.
    """
    DeclaredSymbols = [Symbol('{}'.format(el),commutative=False) for el in symbs]
    Subs = tuple([(el,0) for el in DeclaredSymbols])
    return SimplifytraceMatM(xpr.subs(Subs))
        

# <codecell>

def runningwiththreshold(rges1,rges2,threshold,t0,tf,ts,Assumptions1,Assumptions2,nbpoints=20.) :
    """ 
    does the running of the first set from t0 to t such that t0< t =< threshold
    set the initial values for the second set at t=threshold
    does the running of set 2 from threshold to tf
    In the first version there is no threshold corrections so that all parameters above the scale have to new 
    or define before.
    Warning the initial values of the new parameters that enter at Threshold should be set at this scale
    """
    #Solve the firts set from t0 to 
    if ( threshold - t0 ) /nbpoints < ts :
        print " changing the step to {} for the first section".format((threshold - t0)/nbpoints)
        tts = (threshold - t0)/nbpoints
    rges1.solve_rges(t0,threshold,tts,Assumptions1)
    #the last value for t calculated is not necessarily the same one as the threshold I defined
    tthreshold = rges1.Sol['t'][-1]
    #copy all the parameters in the rges1 into the rges2 that are equal if not found exit and raise an error
    for rg,val in rges1.Sol.items() :
        if not(rg in rges2.labels) and rg != 't':
            print"error {} not in rges2, Check the labels".format(rg)
            exit('error {} not in rges2, Check the labels'.format(rg))
        else :
            #get the position of rg in the label list of rges2
            Idx = [iel for iel,el in enumerate(rges2.labels) if el == rg]
            #copy the value
            if Idx == [] and rg != 't' :
                print "error label {} is not in rges2".format(rg)
                exit('error label {} is not in rges2'.format(rg))
            if rg != 't' :
                rges2.Y0[Idx[0]] = val[-1]
    #Let's adapt the step size again 
    if (tf - tthreshold ) /nbpoints < ts :
        tts = (tf-tthreshold)/ nbpoints
        print " changing the step to {} for the second section".format((tf - tthreshold)/nbpoints)
    #we can now solve the second set of RGEs
    rges2.solve_rges(tthreshold,tf,tts,Assumptions2)       
    #Combine the results
    CombineSol = {}
    for key,val in rges1.Sol.items() :
        CombineSol[key] = np.zeros(len(val)+len(rges2.Sol[key]))
        CombineSol[key][:len(val)] = val
        CombineSol[key][len(val):] = rges2.Sol[key]
    return CombineSol

# <codecell>

def PlotComparison(beta1,beta2,Inputs,moreinfo = '',modelnames=['beta1','beta2']):
    Threshold,t0,tf,ts,Assumptions1,Assumptions2,toplotcombine,toplotSinglet = Inputs
    CombineSol = runningwiththreshold(beta1,beta2,Threshold,t0,tf,ts,Assumptions1,Assumptions2)
    #solve again the SM
    beta1.solve_rges(t0,tf,ts,Assumptions)
    #Plots
    figsize(11,9)
    for ik,k in enumerate(toplotcombine) :
        sx = subplot(len(toplotcombine)+len(toplotSinglet)/2,2,ik+1)
        plt.plot(CombineSol['t'],CombineSol[k],label='beta {} , {} '.format(k,modelnames[0]))
        plt.plot(beta1.Sol['t'],beta1.Sol[k],label='beta {}, {}'.format(k,modelnames[1]))
        #plt.autoscale(tight=True)
        plt.xlim([2,20])
        leg = plt.legend()
        leg.get_frame().set_alpha(0.4)
    for ik,k in enumerate(toplotSinglet) :
        sx = subplot(len(toplotcombine)+len(toplotSinglet)/2,2,ik+1+len(toplotcombine))
        plt.plot(beta2.Sol['t'],beta2.Sol[k],label='beta {} , {}'.format(k,modelnames[1]))
        plt.xlim([2,20])        
        #plt.autoscale(tight=True)
        leg = plt.legend()
        leg.get_frame().set_alpha(0.4)

    plt.suptitle('Different parameters running, {}'.format(moreinfo),y=1.02,fontsize = 14)
    plt.tight_layout()


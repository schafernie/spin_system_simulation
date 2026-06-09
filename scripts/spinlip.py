import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import chi2

def calc_stats(vec):
    n = len(vec)
    if n == 1:
        return np.mean(vec), 0
    return np.mean(vec), np.sqrt(np.var(vec) / (n - 1))



def fitting_routine(fit_function, xdat, ydat, p0=None, yerr=None, printing=True):
    if printing:
        print("***fitting***")
    n = len(xdat)
    if not type(yerr) in [list,np.ndarray]:
        yerr = np.ones(n)
        print("error is set to ones.")
    p, covm = curve_fit(fit_function,xdat,ydat,p0,sigma=yerr) 
    chisq = np.sum(((fit_function(xdat, *p) - ydat)/yerr)**2)
    ndf = n - len(p)
    quality = 1. - chi2.cdf(chisq, ndf)
    chisq = chisq / ndf
    perr = np.sqrt(np.diag(covm)/chisq)      
    if printing:
        print("number of data points: %d" %(n))
        print("dregrees of freedom ndf=%d" %(ndf))
        print("x_min=%f , x_max=%f" %(np.min(xdat),np.max(xdat)))
        print("y_min=%f , y_max=%f" %(np.min(ydat),np.max(ydat)))
        for c in range(0,len(p)):
            print("parameter_%d=%f (+- %f )" %(c,p[c],perr[c]))
        print("q=%f" %(quality)) 
        print("chisq/ndf=%f" %(chisq)) 
    return p, perr, quality
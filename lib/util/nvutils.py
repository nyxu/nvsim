from qutip import *
from datetime import datetime
from time import strftime, localtime



"""get which index of quantum number m in spin number s """
def get_m_index(s,ms):
    #print s,m
    i = s
    n = 0
    while i >= -s :
        if ms == i:
            return int(n)
        i = i - 1
        n = n + 1
    return NaN

def get_ms_set(s):
    mss = []
    i = s
    while i >= -s :
        mss.append(i)
        i = i - 1
        
    return mss

def nvrecord(obj,title = '',dir = '.'):
    time = strftime("%Y%m%d%H%M%S", localtime())
    dir = dir.strip()
    if dir[-1]!='/':
        dir = dir + '/'
    if title == '':
        title = obj.__class__.__name__
    filename = dir + title +time
    qsave(obj,filename)
    print 'object ' + title + 'is saved in ' + filename
    

"""cal mesovle function by a list of argurments,which is for parfor running"""
def mesolve_arglist(arg_array):    
    [Ht, init_dm, evolist, c_ops,op_ops , evo_args,ode_ops,claim] = arg_array
    print claim
    return mesolve(Ht, init_dm, evolist, c_ops,op_ops , evo_args,ode_ops)

def saveHam2Mat(ham,matfile):
    scipy.io.savemat(matfile,{'Hini':ham.data.toarray()})
    
def max_val(mx):
    import sys
    from numbers import Number
    max = -sys.float_info.max
    if isinstance(mx,Number):
        return mx
    for m in mx:
        for v in m:
            if v > max:
                max = v
    return max

"""get the transition gap(frequency) from eigenstate m0 to m1 """
def get_transition_freq(ham,state0,state1):
    gap = -((state0.dag()*ham*state0).tr()-(state1.dag()*ham*state1).tr())
    return gap   

def cal_state_fidelity(s1,s2):
    return abs((s1.dag()*s2).tr())  

from numpy import sqrt
def cal_eigenstate(ham,state,state_discrim_threshold=sqrt(0.5)):
    evals,estates=ham.eigenstates()
    bas = state 
    for s in range(len(estates)):    
        if cal_state_fidelity(bas, estates[s]) > state_discrim_threshold:               
            return estates[s]
         
#     print 'cannot find estates for ',m
    return None

# def cal_fidelity(U1,U2,Ps= None):
#     U1 = U1 / cal_norm(U1)
#     U2 = U2 / cal_norm(U2)
#     if Ps != None:
#         Ua = 0
#         Ub = 0
#         norm =0
#         for P in Ps:
#             Ua = Ua + P*U1*P
#             Ub = Ub + P*U2*P
#             norm = norm + P.tr()
#         return abs((Ua.dag()*Ub).tr())/float(norm)
#     else:
#         return abs((U1.dag()*U2).tr())/float(U1.shape[0])
def cal_fidelity(U1,U2):
#     U1 = U1 / cal_norm(U1)
#     U2 = U2 / cal_norm(U2)

    return abs((U1.dag()*U2).tr())/float(U1.shape[0])

def cal_norm(U):
    return abs((U.dag()*U).tr())/float(U.shape[0])        
            
if __name__=='__main__':
    obj = Qobj()
    nvrecord(obj)
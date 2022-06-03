# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 12:21:05 2020

@author: rmish
"""
import numpy as np 
import scipy.stats as stats
from sklearn.utils import resample
import math
import smtplib
import scipy.stats
from pymoo.model.problem import Problem
from pymoo.algorithms.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.factory import get_termination

def ternary(n,nreq,intvec):
    n = int(n)
    if n == 0:
        return '0'
    nums = []
    while n:
        n, r = divmod(n, 3)
        nums.append(str(r))
    numjoined= ''.join(reversed(nums))
    if intvec == 0:
       while len(numjoined)<(nreq**2):
           numjoined = '0'+numjoined
       tm1 = [int(x)-1 for x in numjoined]
       ternmat = np.reshape(tm1,(nreq,nreq))
    elif intvec == 1:
        ternmat = [int(x)-1 for x in numjoined]    
    return ternmat

def LHS_Sampler(a0,b0, nsamples):
    dim = len(a0) # Number of continuous random variables
    # nsamples = 10
    u = np.zeros(shape=(nsamples, dim))
    samples = np.zeros_like(u)
        
    cut = np.linspace(0, 1, nsamples + 1)
    a = cut[:nsamples]
    b = cut[1:nsamples + 1]
     
    for i in range(dim):
        u[:, i] = stats.uniform.rvs(size=nsamples)
        samples[:, i] = u[:, i] * (b - a) + a
    
    lhs_samples = np.zeros_like(samples)
    for j in range(samples.shape[1]):
        order = np.random.permutation(nsamples)
        lhs_samples[:, j] = samples[order, j]
       
    params = np.zeros_like(lhs_samples)
    for j in range(samples.shape[1]):
        params[:, j] =  lhs_samples[:, j] * (b0[j] - a0[j]) + a0[j]
    
    return params

def Response_Plotter_V1(GM,samp,CAs,lk=[],G_int_vec=[]):
    if lk == []:
        lk = 0
    Os = GM.ortho_nodes
    Is = GM.ind_nodes
    Gconc = samp[:Is]
    dAconc = samp[Is:Is+Os].tolist()    
    dBconc = samp[(Is+Os):(Is+(2*Os))].tolist()
    if G_int_vec == []:        
        binG = bin(int(samp[-1]))[2:]
        while len(binG)<(Os-2):
            binG = '0'+binG
        binGarr = np.array([int(x) for x in binG])
        binGarr[binGarr==0]=-1
        G_int_vec = np.zeros((Is))
        acts = np.array(GM.act_vec)
        for j in range(Os):
            k = j+1
            nn = np.where(acts==k)[0]
            if k == 1 or k == Os: #setting first and last nodes to blocked
                G_int_vec[nn]=-1
            else:
                G_int_vec[nn]=binGarr[k-2] #Setting all other nodes to random values       
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds 
    fstring = 'GdA'+str(Is)
    inC = [0]*Os
    CA = 0
    responseList = []
    for CA in CAs:
        inC[0]=CA
        GM.initial_conditions(dAconc,Gconc,G_int_vec=G_int_vec,dB_added = dBconc,rCin=inC)
        GM.simulate(t_vec1,1,leak=lk)
        G3 = GM.output_concentration[fstring]
        responseList.append(G3)
    return t_vec1, responseList
    

def DR_Plotter_V1(GM,samp,lk=[],maxiter=[],G_int_vec=[],threshold=[]):
    if lk == []:
        lk = 0
    if maxiter == []:
        maxiter = 400
    if threshold == []:
        threshold = 0.98
    Os = GM.ortho_nodes
    Is = GM.ind_nodes
    Gconc = samp[:Is]
    dAconc = samp[Is:Is+Os].tolist()    
    dBconc = samp[(Is+Os):(Is+(2*Os))].tolist()
    if G_int_vec == []:        
        binG = bin(int(samp[-1]))[2:]
        while len(binG)<(Os-2):
            binG = '0'+binG
        binGarr = np.array([int(x) for x in binG])
        binGarr[binGarr==0]=-1
        G_int_vec = np.zeros((Is))
        acts = np.array(GM.act_vec)
        for j in range(Os):
            k = j+1
            nn = np.where(acts==k)[0]
            if k == 1 or k == Os: #setting first and last nodes to blocked
                G_int_vec[nn]=-1
            else:
                G_int_vec[nn]=binGarr[k-2] #Setting all other nodes to random values       
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds 
    fstring = 'GdA'+str(Is)
    GfC = Gconc[-1] #G final max concentration
    inC = [0]*Os
    thresh = 0 #threshold gain (at 90%)
    CA = 0
    resp = np.zeros((maxiter,1))
    dose = np.zeros((maxiter,1))
    counts= 0
    while thresh < threshold and counts < maxiter:
        inC[0]=CA
        GM.initial_conditions(dAconc,Gconc,G_int_vec=G_int_vec,dB_added = dBconc,rCin=inC)
        GM.simulate(t_vec1,1,leak=lk)
        G3 = GM.output_concentration[fstring]
        thresh = G3[-1]/GfC
        # print(counts)
        resp[counts]=G3[-1]/GfC
        dose[counts]=CA/GfC            
        CA += 0.1
        counts+=1
    
    resp = np.trim_zeros(resp)
    dose = dose[:len(resp)] 
    
    return dose,resp,GfC

def Gain_Scorer_V4(GM,samp,lk=[],G_int_vec=[]):
    if lk == []:
        lk = 0
    Os = GM.ortho_nodes
    Is = GM.ind_nodes
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds  
    fstring = 'GdA'+str(Is)
    Gconc = samp[:Is]
    dAconc = samp[Is:Is+Os].tolist()
    dBconc = samp[(Is+Os):(Is+(2*Os))].tolist()
    if G_int_vec == []:            
        binG = bin(int(samp[-1]))[2:]
        while len(binG)<(Os-2):
            binG = '0'+binG
        binGarr = np.array([int(x) for x in binG])
        binGarr[binGarr==0]=-1
        G_int_vec = np.zeros((Is))
        acts = np.array(GM.act_vec)
        for j in range(Os):
            k = j+1
            nn = np.where(acts==k)[0]
            if k == 1 or k == Os: #setting first and last nodes to blocked
                G_int_vec[nn]=-1
            else:
                G_int_vec[nn]=binGarr[k-2] #Setting all other nodes to given samp values   
    GfC = Gconc[-1] #final genelet max concentration
    thresh = 0 #threshold gain (at 90%)
    inC = [0]*Os
    low = 1
    high = Gconc[-1]*10
    counter = 1
    while low + 0.01 <= high: #if the search doesn't find anything in the low/high bounds, make this not happen
        mid = (high + low) / 2
        inC[0]=mid
        GM.initial_conditions(dAconc,Gconc,G_int_vec=G_int_vec,dB_added = dBconc,rCin=inC)
        GM.simulate(t_vec1,1,leak=lk)
        Gout = GM.output_concentration[fstring]
        thresh = Gout[-1]/GfC
        counter +=1
        # print('count:'+str(counter))
        if thresh < 0.89:
            low = mid
        elif thresh > 0.91:
            high = mid
        else:
            return GfC*0.9/mid
            # break
        thresh = 0
    return (GfC*0.9/mid)

def Leak_Scorer_V3(GM,samp,lk=[],G_int_vec=[]): #with standard genelet node stuff
    '''
    Generalized version for any amount of nodes
    '''
    if lk == []:
        lk = 0
    Os = GM.ortho_nodes
    Is = GM.ind_nodes
    Gconc = samp[:Is]
    dAconc = samp[Is:Is+Os].tolist()
    
    dBconc = samp[(Is+Os):(Is+(2*Os))].tolist()
    if G_int_vec==[]:                
        binG = bin(int(samp[-1]))[2:]
        while len(binG)<(Os-2):
            binG = '0'+binG
        binGarr = np.array([int(x) for x in binG])
    
        binGarr[binGarr==0]=-1
        G_int_vec = np.zeros((Is))
        acts = np.array(GM.act_vec)
    
        for j in range(Os):
            k = j+1
            nn = np.where(acts==k)[0]
            if k == 1 or k == Os: #setting first and last nodes to blocked
                G_int_vec[nn]=-1
            else:
                G_int_vec[nn]=binGarr[k-2] #Setting all other nodes to random values
           
    t_vec1 = np.linspace(0,3,1001)*3600 # seconds 
    fstring = 'GdA'+str(Is)
    GfC = Gconc[-1] #G final max concentration
    inC = [0]*Os
    GM.initial_conditions(dAconc,Gconc,G_int_vec=G_int_vec,dB_added = dBconc,rCin=inC)
    GM.simulate(t_vec1,1,leak=lk)
    Gout = GM.output_concentration[fstring]
    leakF = Gout[-1]/GfC
    # print('Leak: ', leakF)
    return leakF

  
def pareto_plotter_v2(AA,lk,numpop,ngen, bounds=[]):
    if bounds == []:        
        G_max = 50
        dA_max = 500
        dB_max = 500
    else:
        G_max = bounds[0]
        dA_max = bounds[1]
        dB_max = bounds[2]
    Os = AA.ortho_nodes
    Is = AA.ind_nodes
    G_int_max = 0# 2**(Os-2) #CHANGE BACK
    # t_vec1 = np.linspace(0,3,1001)*3600 # seconds
    all_max = np.concatenate([G_max*np.ones(Is),dA_max*np.ones(Os),dB_max*np.ones(Os),[G_int_max]])
    cv = np.size(all_max)
    all_min = np.ones(cv)*5
    all_min[Is:Is+Os] = 100 #make sure to have enough activator
    cv = np.size(all_max)
    class topoptimizer(Problem):
    
        def __init__(self, **kwargs):
            super().__init__(n_var=cv,n_obj=2,n_constr=0,xl = all_min, xu=all_max) 
    
        def _evaluate(self, x, out, *args, **kwargs):
            x[:,-1] = np.floor(x[:,-1])
            Gain = np.zeros((numpop,1))
            leak = np.zeros((numpop,1))
            for i in range(len(x)):
                leak[i] = Leak_Scorer_V3(AA,x[i,:],lk=lk)
                # if leak[i] > 0.2:
                #     leak[i] = 1
                #     Gain[i] = 0
                # else:
                Gain[i] = -1*(Gain_Scorer_V4(AA,x[i,:],lk=lk,G_int_vec=list(np.ones(Is)*-1)))
                # Gain[i],leak[i] = AAM.Amplifier_Scorer_v6(AA,x[i,:],lt)
                # print(str(leak[i]))
                # print(str(Gain[i]))
            out["F"]=np.column_stack([Gain,leak])

    p1 = topoptimizer()
    algorithm = NSGA2(pop_size=numpop)
    termination = get_termination("n_gen", ngen)
    res = minimize(p1,algorithm,termination,save_history=True,verbose=True)
    return res

def random_plotter_v2(AA,lk,samps_per_top, bounds = [50,500,500]):
    G_max = bounds[0]
    dA_max = bounds[1]
    dB_max = bounds[2]
    Os = AA.ortho_nodes
    Is = AA.ind_nodes
    G_int_max = 2**(Os-2) 
    all_max = np.concatenate([G_max*np.ones(Is),dA_max*np.ones(Os),dB_max*np.ones(Os),[G_int_max]])
    cv = np.size(all_max)
    all_min = np.ones(cv)*5
    all_min[Is:Is+Os] = 100 #make sure to have enough activator
    all_samp = LHS_Sampler(all_min,all_max,samps_per_top)
    # all_samp[:,-1]=np.floor(all_samp[:,-1]) #rounding the gintvec matrix
    all_samp[:,-1]=np.zeros(len(all_samp[:,-1])) #rounding the gintvec matrix

    again = np.zeros((samps_per_top,1))
    aleak = np.zeros((samps_per_top,1))
    for i in range(samps_per_top):
        again[i]=Gain_Scorer_V4(AA,all_samp[i,:],lk=lk)
        aleak[i]= Leak_Scorer_V3(AA,all_samp[i,:],lk = lk)
        if i%20 == 0:
            print('samps done: ', i)
    return again,aleak, all_samp


def bootstrap_confidence_interval(data,alpha,function):
    """
    To find the bootstrap confidence interval, give the data output (the leak/gain)
     and the 0.05 alpha, and the mean, and let this do its thing
    """
    # Returns lower and upper confidence interval
    B = 19000
    vals = np.empty([B,1], dtype=float)
    N = len(data)
    for i in range(0,B):
        sample = resample(data, n_samples=N)
        vals[i] = function(sample)
    sorted_vals = np.sort(vals,axis=0)
    lower_val_interval = sorted_vals[math.floor(B*alpha/2)]
    upper_val_interval = sorted_vals[math.ceil(B*(1-alpha/2))]
    return [lower_val_interval,upper_val_interval]


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


def sendemail(subject, message,login='simulationpc2021@gmail.com', password='schulmanlab',
              from_addr='simulationpc2021@gmail.com', to_addr_list='rmisha99@gmail.com', 
              cc_addr_list=[],smtpserver='smtp.gmail.com:587'):
        header  = 'From: %s\n' % from_addr
        header += 'To: %s\n' % ','.join(to_addr_list)
        header += 'Subject: %s\n\n' % subject
        message = header + message

        server = smtplib.SMTP(smtpserver)
        server.starttls()
        server.login(login,password)
        server.sendmail(from_addr, to_addr_list, message)
        server.quit()
        
        

        

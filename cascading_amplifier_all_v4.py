# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:14:21 2021

@author: rmish
"""


import numpy as np 
import GeneralGeneletModel_v2 as GAA
import AutoAmplifierModel_v2 as AAM
import matplotlib.pyplot as plt
import pickle
import time

ngen = 150
numpop = 75

def ParetoSaver(pareto,label):
    resF1 = pareto.F
    # resF1 = resF1[resF1[:,1].argsort()]
    resX1 = pareto.X
    RandC1 = np.concatenate((resF1,resX1),axis=1)
    RCorg1 = RandC1[RandC1[:,1].argsort()]
    RCorg1[:,0] = -RCorg1[:,0]
    pickle.dump(RCorg1,open(label,'wb'))
    return resF1

#%% 1 genelet amplifier:
act_vec =  [1,2]
blk_vec =  [1,2]
prod_vec = [2,1]
indc_vec = [0,0]
G_int_vec = [-1,-1]

AA = GAA.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
AA.plot_topology(show_rnas=0)

t1 = time.time()
pareto1 = AAM.pareto_plotter_v2(AA,0.01,numpop,ngen,bounds=[100,1000,1000])
t2 = time.time()

resF1 = ParetoSaver(pareto1,'1geneletamp_'+str(numpop)+'_'+str(ngen)+'.obj')

#%% 2 genelet amplifier:
act_vec =  [1,2,3]
blk_vec =  [1,2,3]
prod_vec = [2,3,0]
indc_vec = [0,0,0]
G_int_vec = [-1,-1,-1]
AA = GAA.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
AA.plot_topology(show_rnas=0)

t3 = time.time()
pareto2 = AAM.pareto_plotter_v2(AA,0.01,numpop,ngen,bounds=[100,1000,1000])
t4 = time.time()
resF2 = ParetoSaver(pareto2,'2geneletamp_'+str(numpop)+'_'+str(ngen)+'.obj')


#%% 3 genelet amplifier:
act_vec =  [1,2,3,4]
blk_vec =  [1,2,3,4]
prod_vec = [2,3,4,0]
indc_vec = [0,0,0,0]
G_int_vec = [-1,-1,-1,-1]
AA = GAA.GeneletNetwork(act_vec,prod_vec,indc_vec,blk_vec)
AA.plot_topology(show_rnas=0)

t5 = time.time()
pareto3 = AAM.pareto_plotter_v2(AA,0.01,numpop,ngen,bounds=[100,1000,1000])
t6 = time.time()
resF3 = ParetoSaver(pareto3,'3geneletamp_'+str(numpop)+'_'+str(ngen)+'.obj')


#%%Plotting and saving:

print('time1: '+str(t2-t1))
print('time2: '+str(t4-t3))
print('time3: '+str(t6-t5))
    
plt.figure(dpi=150)
plt.scatter(resF1[:,1],-resF1[:,0],color='r')
plt.plot(resF1[:,1],-resF1[:,0],color='r',linewidth=4)
plt.scatter(resF2[:,1],-resF2[:,0],color='b')
plt.plot(resF2[:,1],-resF2[:,0],color='b',linewidth=4)
plt.scatter(resF3[:,1],-resF3[:,0],color='g')
plt.plot(resF3[:,1],-resF3[:,0],color='g',linewidth=4)
ax1 = plt.gca()
ax1.set_xlim(0,0.2)
# ax1.set_ylim(0,50)
fs = 13
plt.title('Cascading Amplifier Pareto Fronts',fontsize=15,weight='bold')
plt.ylabel('Normalized Amplification',fontsize=fs,weight='bold')
plt.xlabel('Normalized Leak',fontsize=fs,weight='bold')
plt.xticks(fontsize=fs,weight='bold')
plt.yticks(fontsize=fs,weight='bold')
plt.locator_params(nbins=6)











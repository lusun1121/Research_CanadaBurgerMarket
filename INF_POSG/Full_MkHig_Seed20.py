# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 09:56:20 2022

@author: lusun8825
"""

#import scipy.io
import numpy as np
import pandas as pd
#import itertools
#import time
from scipy.special import logsumexp
from scipy.special import softmax

mk=3
seed_ini = 20
theta23 = np.array([0.8993, 0.0080]+ [0.9172])#, 0.1,0.9,0.1])#,0.1,0.1])
thetaS = np.array([1.8028, 1.9129, -0.3808 , -0.0252 ,-0.0004, 0.0174, 0.0052 ,9.0096,
                   1.3930,1.4743, -0.3980, -0.0144,0.0083,0.0089, 0.0250, 7.1176])
# thetaS = np.array([1.9160,2.1589 ,-0.5025,-0.0981,0.0511,0.0562,0.0516,12.6968,
#                    0.9989,1.4974,-0.4957,-0.1537,0.0288,0.0586,0.0762,7.5245]) #gdp posg full mk=2

#path = '/scratch/user/lusun8825/rust_hidden/'
file_name = 'G:/My Drive/2022_Fall/Research_Burger/INF_POSG/'


dataFull = pd.read_csv(file_name + 'dataFull.csv')
TimePeriod = 35#len(np.unique(data[:,1]))
SamplePath = int(len(dataFull)/TimePeriod)
MT = len(dataFull)          #% Size of the dataset (M x T = 400 x 35 = 14000).#clusterid = np.int_(data[:,0])      % Unique market ID.

tertile = dataFull['tertile'].to_numpy(dtype=int)      #% Market type initial guess (from 130913_alaToivanenWaterson3_Fixedeffect3quantile.csv/dta)
#% Number of own shops (in state space; capped at 3), from the perspective of each firm
Ni_aw = dataFull['Ni_aw'].to_numpy(dtype=int) 
Ni_bk = dataFull['Ni_bk'].to_numpy(dtype=int) 
Ni_hvy = dataFull['Ni_hvy'].to_numpy(dtype=int) 
Ni_mcd = dataFull['Ni_mcd'].to_numpy(dtype=int) 
Ni_wdy = dataFull['Ni_wdy'].to_numpy(dtype=int) 
#% Number of rival shops (in state space; capped at 3), from the perspective of each firm
Nj_aw =dataFull['Nj_aw'].to_numpy(dtype=int)  
Nj_bk = dataFull['Nj_bk'].to_numpy(dtype=int) 
Nj_hvy =dataFull['Nj_hvy'].to_numpy(dtype=int) 
Nj_mcd = dataFull['Nj_mcd'].to_numpy(dtype=int) 
Nj_wdy = dataFull['Nj_wdy'].to_numpy(dtype=int) 
#% Number of pop,val,inc (in state space; capped at 3) 
disc_pop = dataFull['disc_pop'].to_numpy(dtype=int) 
disc_val = dataFull['disc_val'].to_numpy(dtype=int) 
disc_inc = dataFull['disc_inc'].to_numpy(dtype=int) 

stateID_aw = Ni_aw*4**4+Nj_aw*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_bk = Ni_bk*4**4+Nj_bk*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_hvy = Ni_hvy*4**4+Nj_hvy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_mcd = Ni_mcd*4**4+Nj_mcd*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_wdy = Ni_wdy*4**4+Nj_wdy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_other = np.hstack([stateID_aw,stateID_bk,stateID_hvy,stateID_wdy])

ai_aw = dataFull['ai_aw'].to_numpy(dtype=int) 
ai_bk =dataFull['ai_bk'].to_numpy(dtype=int) 
ai_hvy = dataFull['ai_hvy'].to_numpy(dtype=int) 
ai_mcd = dataFull['ai_mcd'].to_numpy(dtype=int) 
ai_wdy = dataFull['ai_wdy'].to_numpy(dtype=int) 
ai_other = np.hstack([ai_aw,ai_bk,ai_hvy,ai_wdy])

numstates = 5 # ni,nj,pop,val,inc
states_full = np.arange(4**numstates)
Ni = states_full//(4**(numstates-1))                       #Number of own outlets: {0,1,2,3+}.
Nj = states_full%(4**(numstates-1))//(4**(numstates-2))    #Number of rival outlets: {0,1,2,3+}.
dz1 = states_full%(4**(numstates-2))//(4**(numstates-3))   #Discretized population state: {0,1,2,3}.
dz2 = states_full%(4**(numstates-3))//(4**(numstates-4))   #Discretized property value state: {0,1,2,3}.
dz3 = states_full%(4**(numstates-4))                       #Discretized income state: {0,1,2,3}.
#print(np.sum(np.abs(Ni*4**4+Nj*4**3+dz1*4**2+dz2*4**1+dz3-states_full)))
states = np.stack([Ni,Nj,dz1,dz2,dz3],axis=1)

# disc_inf = dataFull['disc_inf'].to_numpy(dtype=int) 
# disc_gdp = dataFull['disc_gdp'].to_numpy(dtype=int) 
disc_gdp = dataFull['disc_inf'].to_numpy(dtype=int) 

#%%

I = 2 #MCD and Other homogenous
J = 3 #Action: 0,1,2
NS = 4**(numstates)
NX = 15
obser = 2
NZ = obser

F_Ni = np.zeros([I,J,NX,NS,NS])#;
F_Ni[0,0] = np.kron(np.array([[1]*4**4+[0]*4**4 + [0] *4**4+[0]*4**4]*4**4 +\
                          [[1]*4**4+[0]*4**4 + [0] *4**4+[0]*4**4]*4**4 + \
                          [[0]*4**4+[1]*4**4 + [0] *4**4+[0]*4**4]*4**4 + \
                          [[0]*4**4+[0]*4**4 + [1] *4**4+[0]*4**4]*4**4 ),np.ones([NX,1,1]))
F_Ni[0,1] =np.kron(np.array([[1]*4**4+[0]*4**4 + [0] *4**4+[0]*4**4]*4**4 +\
                          [[0]*4**4+[1]*4**4 + [0] *4**4+[0]*4**4]*4**4 + \
                          [[0]*4**4+[0]*4**4 + [1] *4**4+[0]*4**4]*4**4 + \
                          [[0]*4**4+[0]*4**4 + [0] *4**4+[1]*4**4]*4**4 ),np.ones([NX,1,1]))
F_Ni[0,2] = np.kron(np.array([[0]*4**4+[1]*4**4 + [0] *4**4+[0]*4**4]*4**4 +\
                          [[0]*4**4+[0]*4**4 + [1] *4**4+[0]*4**4]*4**4 + \
                          [[0]*4**4+[0]*4**4 + [0] *4**4+[1]*4**4]*4**4 + \
                          [[0]*4**4+[0]*4**4 + [0] *4**4+[1]*4**4]*4**4 ),np.ones([NX,1,1]))
F_Ni[1] = F_Ni[0]

disc_pop_mx = disc_pop[tertile==mk].reshape([-1,TimePeriod])
disc_pop_bk = disc_pop_mx[:,0:TimePeriod-1].reshape(-1)
disc_pop_fw = disc_pop_mx[:,1:TimePeriod].reshape(-1)
disc_val_mx = disc_val[tertile==mk].reshape([-1,TimePeriod])
disc_val_bk = disc_val_mx[:,0:TimePeriod-1].reshape(-1)
disc_val_fw = disc_val_mx[:,1:TimePeriod].reshape(-1)
disc_inc_mx = disc_inc[tertile==mk].reshape([-1,TimePeriod])
disc_inc_bk = disc_inc_mx[:,0:TimePeriod-1].reshape(-1)
disc_inc_fw = disc_inc_mx[:,1:TimePeriod].reshape(-1)

#% F_dz: Make simple & intuitive (4x4) versions for exposition purposes
fz1_4x4 = np.zeros([4,4])#;
fz2_4x4 = np.zeros([4,4])#;
fz3_4x4 = np.zeros([4,4])#;

for x0 in range(4):
    for x1 in range(4):
        numer1 = np.sum((disc_pop_fw == x1) * (disc_pop_bk == x0))#;
        numer2 = np.sum((disc_val_fw == x1) * (disc_val_bk == x0))#;
        numer3 = np.sum((disc_inc_fw == x1) * (disc_inc_bk == x0))#;
        denom1 = np.sum((disc_pop_bk == x0))#;
        denom2 = np.sum((disc_val_bk == x0))#;
        denom3 = np.sum((disc_inc_bk == x0))#;

        fz1_4x4[x0,x1] = numer1 / denom1#;
        fz2_4x4[x0,x1] = numer2 / denom2#;
        fz3_4x4[x0,x1] = numer3 / denom3#;

# for i in range(4):
#   print("&{}&{:.4f}&{:.4f}&{:.4f}&{:.4f}\\\\".format(i,fz1_4x4[i][0],fz1_4x4[i][1],fz1_4x4[i][2],fz1_4x4[i][3]))
# print('------------------------------')
# for i in range(4):
#   print("&{}&{:.4f}&{:.4f}&{:.4f}&{:.4f}\\\\".format(i,fz2_4x4[i][0],fz2_4x4[i][1],fz2_4x4[i][2],fz2_4x4[i][3]))
# print('------------------------------')
# for i in range(4):
#   print("&{}&{:.4f}&{:.4f}&{:.4f}&{:.4f}\\\\".format(i,fz3_4x4[i][0],fz3_4x4[i][1],fz3_4x4[i][2],fz3_4x4[i][3]))
# print('------------------------------')
# print(fz1_4x4)
# print(fz2_4x4)
# print(fz3_4x4)

F_dz1 = np.kron((fz1_4x4[dz1,:])[:,dz1],np.ones([I,J,NX,1,1]))#; % The transition pattern is common across (i,j,mtype)
F_dz2 = np.kron((fz2_4x4[dz2,:])[:,dz2],np.ones([I,J,NX,1,1]))#; % The transition pattern is common across (i,j,mtype)
F_dz3 = np.kron((fz3_4x4[dz3,:])[:,dz3],np.ones([I,J,NX,1,1]))#; % The transition pattern is common across (i,j,mtype)


dim_hstate = 2
dim_observe = obser
def Dynamic(theta23):
  trans_val = theta23[0:2]
  obser_val = theta23[2:]

  trans = np.zeros([dim_hstate,dim_hstate]) #s_t^h,s_{t+1}^h
  trans[:,0] = trans_val
  trans[:,1] = 1-trans[:,0]
  obser = np.zeros([dim_hstate,dim_observe])
  obser[:,0] = [obser_val,1-obser_val]
  #obser[:,0] = obser_val
  obser[:,1] = 1-obser[:,0]
  return trans,obser
#print(Dynamic(np.array([0.95, 0.05]+ [0.125, 0.125])))#,0.125]+[0.1])))#,0.1,0.1])))  

def SigmaLambda(theta23,z_old = None,x_old=None,T=None,num_discrete = NX):
  trans,obser = Dynamic(theta23)
  if T==None: #generate function for Q function
    x_old =  np.linspace(0,1,num=num_discrete)
    x_new = np.zeros([num_discrete,dim_observe]) #x_old,z_old,z_new
    sigma = np.zeros([num_discrete,dim_observe])

    for z in range(dim_observe):
      #for z_prime in range(dim_observe):
        x_temp = (x_old*trans[0,0] + (1-x_old) * trans[1,0])*obser[0,z]#,z_prime]
        sigma[:,z] = x_temp + (x_old*trans[0,1] + (1-x_old) * trans[1,1])*obser[1,z]#,z_prime]
        sigma_nonzero = np.where(sigma[:,z]!=0)
        x_new[:,z][sigma_nonzero] = x_temp[sigma_nonzero]/sigma[:,z][sigma_nonzero]
  elif T==1:
    x_temp = (x_old*trans[0,0]+(1-x_old)*trans[1,0])*obser[0,z_old]
    sigma = x_temp + (x_old*trans[0,1]+(1-x_old)*trans[1,1])*obser[1,z_old]
    #sigma_nonzero = np.where(sigma!=0)
    x_new = x_temp/sigma#[sigma_nonzero]/sigma[sigma_nonzero]
  else: #generate whole blief in recover process
    x_new = np.zeros([T,len(x_old)])
    sigma = np.zeros([T-1,len(x_old)])
    x_new[0] = x_old
    for t in range(T-1):
      #z = z_old[t,:]
      z_prime = z_old[t+1,:]
      x_temp = (x_new[t]*trans[0,0]+(1-x_new[t])*trans[1,0])*obser[0,z_prime]
      sigma[t] = x_temp + (x_new[t]*trans[0,1]+(1-x_new[t])*trans[1,1])*obser[1,z_prime]
      sigma_nonzero = np.where(sigma[t]!=0)
      x_new[t+1][sigma_nonzero] = x_temp[sigma_nonzero]/sigma[t][sigma_nonzero]

  belief_f = np.floor(x_new*(num_discrete-1))/(num_discrete-1)
  belief_c = np.ceil(x_new*(num_discrete-1))/(num_discrete-1)
  iterpolate = np.zeros(belief_f.shape)
  iterpolate[np.where((belief_f-belief_c)!=0)] = (x_new-belief_c)[np.where((belief_f-belief_c)!=0)]/(belief_f-belief_c)[np.where((belief_f-belief_c)!=0)]      

  return sigma,x_new,[iterpolate,np.int_(belief_f*(num_discrete-1)),np.int_(belief_c*(num_discrete-1))]

# theta23 = np.array([0.95, 0.05,0.125, 0.125])#,0.125]+[0.1])#,0.1,0.1])
# sigma,_,(iterp,iterf,iterc) = SigmaLambda(theta23)

F_fix = F_Ni * F_dz1 * F_dz2 *F_dz3
# % Calculate F from P (& global transition matrices F_Ni, F_dz1, F_dz2) 
def updateF(P):
  #% For a given combination of (Ni, dz1, dz2, mtype), make the state indexes (x, x0, x1, x2, x3)
  #% Beliefs: mapping from today's state (Ni,Nj,dz1,dz2) to tomorrow's Nj = {0,1,2,3}
  fnj_mcd = np.zeros([NX,NS,4])#;      % McDonald's belief
  fnj_other = np.zeros([NX,NS,4])#;    % Other 4 chains beliefs

  nj = 0#;         % If today's # of rivals = 0
  x =  (4**4)*Ni + (4**3)*nj + (4**2)*dz1 + 4*dz2+dz3#;  % index of today's own state
  x0 = (4**4)*0 + (4**3)*Ni + (4**2)*dz1 + 4*dz2+dz3#;  % index of today's state for rivals with 0 shops

  fnj_mcd[:,x,0] = P[1,1][:,x0]**4#;                               % Prob(Nj'= 0) from McD's perspective
  fnj_mcd[:,x,1] = 4 * P[1,2][:,x0] * (P[1,1][:,x0]**3)#;       % Prob(Nj'= 1) from McD's perspective
  fnj_mcd[:,x,2] = 6 * (P[1,2][:,x0]**2) * (P[1,1][:,x0]**2)#;   % Prob(Nj'= 2) from McD's perspective
  fnj_mcd[:,x,3] = 4 * (P[1,2][:,x0]**3) * P[1,1][:,x0]#;       % Prob(Nj'= 3) from McD's perspective
  
  fnj_other[:,x,0] = P[0,1][:,x0] * (P[1,1][:,x0]**3)#;         % Prob(Nj'= 0) etc. from Others' perspectives
  fnj_other[:,x,1] = P[0,2][:,x0] * (P[1,1][:,x0]**3) + P[0,1][:,x0] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2)#;  
  fnj_other[:,x,2] = P[0,2][:,x0] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2) + P[0,1][:,x0] * 3 * (P[1,2][:,x0]**2) * P[1,1][:,x0]
  fnj_other[:,x,3] = P[0,2][:,x0] * 3 * (P[1,2][:,x0]**2) * P[1,1][:,x0] + P[0,1][:,x0] * (P[1,2][:,x0]**3)
  
  nj = 1#;         % If today's # of rivals = 1                               
  x =  (4**4)*Ni + (4**3)*nj + (4**2)*dz1 +4*dz2+dz3#;  % index of today's own state
  x0 = (4**4)*0 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 0 shops
  x1 =  (4**4)*1 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 1 shops
  
  fnj_mcd[:,x,0] = P[1,0][:,x1] * (P[1,1][:,x0]**3)#;
  fnj_mcd[:,x,1] = P[1,1][:,x1] * (P[1,1][:,x0]**3)#;
  fnj_mcd[:,x,2] = P[1,2][:,x1] * (P[1,1][:,x0]**3) + P[1,1][:,x1] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2)#;
  fnj_mcd[:,x,3] = P[1,2][:,x1] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2) + P[1,1][:,x1] * 3 * (P[1,2][:,x0]**2) * P[1,1][:,x0]
  
  fnj_other[:,x,0] = .5 * (P[0,0][:,x1] * P[1,1][:,x0]**3) + .5 * (P[0,1][:,x0] * P[1,0][:,x1] * (P[1,1][:,x0]**2))
  fnj_other[:,x,1] = .5 * (P[0,1][:,x1] * P[1,1][:,x0]**3) + .5 * (P[0,1][:,x0] * P[1,1][:,x1] * (P[1,1][:,x0]**2))
  fnj_other[:,x,2] = .5 * (P[0,2][:,x1] * (P[1,1][:,x0]**3) + P[0,1][:,x1] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2))+\
                              .5 * (P[0,2][:,x0] * P[1,1][:,x1] * (P[1,1][:,x0]**2) +\
                        P[0,1][:,x0] * (P[1,2][:,x1] * (P[1,1][:,x0]**2) + P[1,1][:,x1] * 2 * P[1,2][:,x0] * P[1,1][:,x0]))
  fnj_other[:,x,3] = .5 * (P[0,2][:,x1] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2)\
                        + P[0,1][:,x1] * 3 * (P[1,2][:,x0]**2) * P[1,1][:,x0])\
                        + .5 * (P[0,2][:,x0] * P[1,2][:,x1] * (P[1,1][:,x0]**2)\
                      + P[0,2][:,x0] * P[1,1][:,x1] * 2 * P[1,2][:,x0] * P[1,1][:,x0]\
                      + P[0,1][:,x0] * P[1,2][:,x1] * 2 * P[1,2][:,x0] * P[1,1][:,x0]\
                      + P[0,1][:,x0] * P[1,1][:,x1] * (P[1,2][:,x0]**2))

  nj = 2#;         % If today's # of rivals = 2
  x = (4**4)*Ni + (4**3)*nj + (4**2)*dz1 +4*dz2+dz3#;  % index of today's own state
  x0 = (4**4)*0 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 0 shops
  x1 = (4**4)*1 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 1 shops
  x2 = (4**4)*2 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 2 shops
  
  fnj_mcd[:,x,0] = (P[1,0][:,x1]**2) * (P[1,1][:,x0]**2)
  fnj_mcd[:,x,1] = (1/3) * (P[1,0][:,x2] * (P[1,1][:,x0]**3)) + (2/3) * (P[1,0][:,x1] * P[1,1][:,x1] * (P[1,1][:,x0]**2))
  fnj_mcd[:,x,2] = (1/3) * (P[1,1][:,x2] * (P[1,1][:,x0]**3)) + (2/3) * ((P[1,1][:,x1]**2) * (P[1,1][:,x0]**2))
  fnj_mcd[:,x,3] = (1/3) * (P[1,2][:,x2] * (P[1,1][:,x0]**3)  + (P[1,1][:,x2] * 3 * P[1,2][:,x0] * (P[1,1][:,x0]**2)))
  
  fnj_other[:,x,0] = .5 * 0 + .5 * (P[0,0][:,x1] * P[1,0][:,x1] * (P[1,1][:,x0]**2))
  fnj_other[:,x,1] = .5 * (P[0,0][:,x2] * (P[1,1][:,x0]**3))\
      + .5 * (P[0,0][:,x1] * P[1,1][:,x1] * (P[1,1][:,x0]**2)\
      + P[0,1][:,x1] * P[1,0][:,x1] * (P[1,1][:,x0]**2))
  fnj_other[:,x,2] = .5 * (P[0,1][:,x2] * (P[1,1][:,x0]**3))\
      + .5 * (P[0,1][:,x1] * P[1,1][:,x1] * (P[1,1][:,x0]**2))
  fnj_other[:,x,3] = .5 * (P[0,2][:,x2] * (P[1,1][:,x0]**3)\
      + P[0,1][:,x2] * 3 * P[1,2][:,x0] * P[1,1][:,x0])\
      + .5 * (P[0,2][:,x1] * P[1,1][:,x1] * (P[1,1][:,x0]**2)\
      + P[0,1][:,x1] * P[1,2][:,x1] * (P[1,1][:,x0]**2)\
      + P[0,1][:,x1] * P[1,1][:,x1] * 2 * P[1,2][:,x0] * P[1,1][:,x0])
  
  nj = 3#;         % If today's # of rivals = 3
  x = (4**4)*Ni + (4**3)*nj + (4**2)*dz1 +4*dz2+dz3#;  % index of today's own state
  x0 = (4**4)*0 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 0 shops
  x1 = (4**4)*1 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 1 shops
  x2 = (4**4)*2 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 2 shops
  x3 = (4**4)*3 + (4**3)*Ni + (4**2)*dz1 +4*dz2+dz3#;  % index of today's state for rivals with 3 shops
  
  fnj_mcd[:,x,0] = (2/3) * (P[1,0][:,x1]**3) * P[1,1][:,x0]
  fnj_mcd[:,x,1] = (2/3) * (P[1,0][:,x1]**2) * P[1,1][:,x1] + (1/3) * P[1,0][:,x2] * P[1,0][:,x1] * (P[1,1][:,x0]**2)
  fnj_mcd[:,x,2] = (2/3) * P[1,0][:,x1] * (P[1,1][:,x1]**2) + (1/3) * (P[1,0][:,x2] * P[1,1][:,x1] * (P[1,1][:,x0]**2)\
      + P[1,1][:,x2] * P[1,0][:,x1] * (P[1,1][:,x0]**2))
  fnj_mcd[:,x,3] = (2/3) * (P[1,1][:,x1]**3) * P[1,1][:,x0] + (1/3) * P[1,1][:,x2] * P[1,1][:,x1] * (P[1,1][:,x0]**2)
  
  fnj_other[:,x,0] = .25 * (P[0,0][:,x1] * (P[1,0][:,x1]**2) * P[1,1][:,x0])
  fnj_other[:,x,1] = .25 * 0 + .25 * (P[0,0][:,x2] * P[1,0][:,x1] * (P[1,1][:,x0]**2))\
      + .25 * (P[0,0][:,x1] * 2 * P[1,0][:,x1] * P[1,1][:,x1] * P[1,1][:,x0]\
      + P[0,1][:,x1] * (P[1,0][:,x1]**2) * P[1,1][:,x0])\
      + .25 * (P[0,0][:,x1] * P[1,0][:,x2] * (P[1,1][:,x0]**2))
  fnj_other[:,x,2] = .25 * (P[0,0][:,x3] * (P[1,1][:,x0]**3))\
      + .25 * (P[0,0][:,x2] * P[1,1][:,x1] * (P[1,1][:,x0]**2)\
      + P[0,1][:,x2] * P[1,0][:,x1] * (P[1,1][:,x0]**2))\
      + .25 * (P[0,0][:,x1] * (P[1,1][:,x1]**2) * P[1,1][:,x0]\
      + P[0,1][:,x1] * 2 * P[1,0][:,x1] * P[1,1][:,x1] * P[1,1][:,x0])\
      + .25 * (P[0,0][:,x1] * P[1,1][:,x2] * (P[1,1][:,x0]**2)\
      + P[0,1][:,x1] * P[1,0][:,x2] * (P[1,1][:,x0]**2))
  fnj_other[:,x,3] = .25 * (P[0,1][:,x3] * (P[1,1][:,x0]**3))\
      + .25 * (P[0,1][:,x2] * P[1,1][:,x1] * (P[1,1][:,x0]**2))\
      + .25 * (P[0,1][:,x1] * (P[1,1][:,x1]**2) * P[1,1][:,x0])\
      + .25 * (P[0,1][:,x1] * P[1,1][:,x2] * (P[1,1][:,x0]**2))

  F_Nj = np.ones([I,J,NX,NS,NS])
  Denom_mcd = np.matmul(fnj_mcd,np.ones([NX,4,4]))#np.kron(np.reshape(np.sum(fnj_mcd,axis=2),[NX,-1,1]),np.ones([1,1,4]))#;         % Sum of each row
  fnj_mcd = fnj_mcd / Denom_mcd#;                 % Make sure each row sums up to 1
  fnj_mcd[np.isnan(fnj_mcd)] = 0#;                    % Replace NaN with 0
  F_Nj[0,0] = fnj_mcd[:,:,Nj]

  Denom_other = np.matmul(fnj_other,np.ones([NX,4,4]))#np.kron(np.reshape(np.sum(fnj_other,axis=2),[NX,-1,1]),np.ones([1,1,4]))#;             % Sum of each row
  fnj_other = fnj_other / Denom_other#;                   % Make sure each row sums up to 1
  fnj_other[np.isnan(fnj_other)] = 0#;                        % Replace NaN with 0
  F_Nj[1,0] = fnj_other[:,:,Nj]
  #F_Nj[2:5,j,:,:] = np.stack([F_Nj[1,j,:,:] for copy in range(3)],axis=0)#;    % Other 4 chains are symmetric

  for j in range(1,J):    
      F_Nj[:,j] = F_Nj[:,0]    
  return F_Nj*F_fix


def updatePi(Q):
  P = np.zeros([I,J,NX,NS])#; player, action, states, Market Type
  P[0] = softmax(Q[0],axis=0)
  P[1] = softmax(Q[1],axis=0)
  P[np.isnan(P)] = 0#;
  return P

def updateQ(thetaS,Q,F,sigmaP,iterpP,iterf,iterc,Q_new,beta = 0.9):

  gamma = 0.5772
    
  v_mcd = np.zeros([J,NX,obser,NS])
  v_other = np.zeros([J,NX,obser,NS])

  v_mcd= iterpP*Q[0][:,iterf,:] + (1-iterpP)*Q[0][:,iterc,:] # Agent, action,NX_new(NX_old,NZ_old,NZ_new), NZ_new,NS_new
  v_other = iterpP*Q[1][:,iterf,:] + (1-iterpP)*Q[1][:,iterc,:]
  #v_mcd = np.stack([iterpP*Q[0,j][iterf]+(1-iterpP)*Q[0,j][iterc] for j in range(J)],axis=0)
  #v_other = np.stack([iterpP*Q[1,j][iterf]+(1-iterpP)*Q[1,j][iterc] for j in range(J)],axis=0)
  v_mcd = gamma + logsumexp(v_mcd,axis=0) #x_old,z_old,z_new,NS_new
  v_other = gamma + logsumexp(v_other,axis=0) #x_old,z_old,z_new,NS_new
  v_mcd = np.sum(sigmaP*v_mcd,axis=1).reshape([1,NX,NS,1]) #x_old,z_old,NS_new
  v_other = np.sum(sigmaP*v_other,axis=1).reshape([1,NX,NS,1])
  v_mcd = np.kron(v_mcd,np.ones([J,1,1,1])) #a_old,x_old,NZ_old,NS_new,1
  v_other = np.kron(v_other,np.ones([J,1,1,1])) #a_old,x_old,NZ_old,NS_new,1

  v_temp_mcd = np.matmul(F[0],v_mcd)[:,:,:,0] # action_old,x_old,nz_old, ns_old (ns_new \times ns_new)
  v_temp_other = np.matmul(F[1],v_other)[:,:,:,0] 

  Q_new[0] = Q_new[0] + beta*(v_temp_mcd) 
  Q_new[1] = Q_new[1] + beta*(v_temp_other)
  return Q_new


def ValueIteration(thetaS,sigmaP,iterpP,iterf,iterc,error = 1e-2,beta = 0.9):
  Pi0 = 1/3*np.ones([I,J,NX,NS])
  F0= updateF(Pi0)
  Q0 = np.zeros([I,J,NX,NS])
  #Pi_old = Pi0
  F_old = F0
  Q_old = Q0

  x_old = np.linspace(0,1,num=NX).reshape([-1,1])
  R_new = np.zeros([I,J,NX,NS])
  thetaS_mcd = np.reshape(np.array([thetaS[2],thetaS[3],thetaS[4],thetaS[5],thetaS[6],thetaS[7]]),[-1,1])
  thetaS_other = np.reshape(np.array([thetaS[10],thetaS[11],thetaS[12],thetaS[13],thetaS[14],thetaS[15]]),[-1,1])
  Ni_temp = (Ni).reshape([1,-1])

  r_mcd_a2 = np.stack([Ni *Ni,Ni*Nj, Ni*dz1,Ni*dz2,Ni*dz3,-1* np.ones(NS)],axis=1).dot(thetaS_mcd)
  r_mcd_a1 = np.stack([Ni *Ni,Ni*Nj, Ni*dz1,Ni*dz2,Ni*dz3,  0* np.ones(NS)],axis=1).dot(thetaS_mcd)

  R_new[0,2] =   np.ones([NX,1]).dot(r_mcd_a2.T)+\
                thetaS[0]*    x_old.dot(Ni_temp)+\
                thetaS[1]*(1-x_old).dot(Ni_temp)
  R_new[0,1] =   np.ones([NX,1]).dot(r_mcd_a1.T)+\
                thetaS[0]*    x_old.dot(Ni_temp)+\
                thetaS[1]*(1-x_old).dot(Ni_temp)
  R_new[0,0] = R_new[0,1,:,:].copy()

  r_other_a2 = np.stack([Ni *Ni,Ni*Nj,Ni*dz1,Ni*dz2,Ni*dz3,  -1* np.ones(NS)],axis=1).dot(thetaS_other)
  r_other_a1 = np.stack([Ni *Ni,Ni*Nj,Ni*dz1,Ni*dz2,Ni*dz3,  0* np.ones(NS)],axis=1).dot(thetaS_other)

  R_new[1,2] =  np.ones([NX,1]).dot(r_other_a2.T)+\
                thetaS[8]*    x_old.dot(Ni_temp)+\
                thetaS[9]*(1-x_old).dot(Ni_temp)
  R_new[1,1] =  np.ones([NX,1]).dot(r_other_a1.T)+\
                thetaS[8]*    x_old.dot(Ni_temp)+\
                thetaS[9]*(1-x_old).dot(Ni_temp)
  R_new[1,0] = R_new[1,1].copy()
  error1 = np.inf
  for iterl in range(1000):
    #time_start = time.time()
    Q_new = updateQ(thetaS,Q_old,F_old,sigmaP,iterpP,iterf,iterc,Q_new = R_new.copy(),beta=beta) #I, J, NS
    Pi_new = updatePi(Q_new)#I,J,NS
    #time_start1 = time.time()
    #print(time_start1-time_start)
    F_new = updateF(Pi_new) # I, J, NS,NS
    #print(time.time()-time_start1)
    if iterl>0:
      error1 = np.max(np.abs(Q_new-Q_old))/np.max(np.abs(Q_old))
    #error2 = np.max(np.abs(Pi_new-Pi_old))/np.max(np.abs(Pi_old))
    #error3 = np.max(np.abs(F_new-F_old))/np.max(np.abs(F_old))
    #print('iteration {} --> e1:{:.4f}({:.4f}), e2:{:.4f}({:.4f}), e3:{:.4f}({:.4f}), time:{}'.format(iter,error1,np.max(np.abs(Q_old)),error2,np.max(np.abs(Pi_old)),error3,np.max(np.abs(F_old)),time.time()-time_start))
    #print('iteration {} --> e1:{},time:{}'.format(iter,error1,time.time()-time_start))

    if error1<error :#and error2<error and error3<error:
      #print('iteration {} --> e1:{},e2:{},e3:{}'.format(iter,error1,error2,error3))
      break
    Q_old = Q_new.copy()
    #Pi_old = Pi_new.copy()
    F_old = F_new.copy()
  return Q_new, Pi_new,F_new


# theta23 = np.array([0.8993, 0.0080]+ [0.9172])#, 0.1,0.9,0.1])#,0.1,0.1])
# sigma,_,(iterp,iterf,iterc) = SigmaLambda(theta23)
# sigmaP = np.kron(sigma.reshape([NX,obser,1]),np.ones([1,1,NS]))
# iterpP = np.kron(iterp.reshape([1,NX,obser,1]),np.ones([J,1,1,NS]))

# thetaS = np.array([2.4672,3.9488,-2.5972,-0.4203,0.0931,0.1513,0.0730,12.5872,
#                    1.1637,1.3230,-0.54,-0.1731,0.0596,0.0413,0.0116,7.6820])
# # #thetaS = np.zeros(5*2)
# Q_new, Pi_new,F_new = ValueIteration(thetaS,sigmaP,iterpP,iterf,iterc)
#%%
tertile_other = np.hstack([tertile,tertile,tertile,tertile])

disc_gdp_other = np.hstack([disc_gdp for i in range(4)])
disc_gdp_mcd = disc_gdp.copy()
z_data = disc_gdp[tertile==mk].reshape([-1,TimePeriod]).T
x_data = (1/obser)*np.ones(z_data.shape[1])

from scipy.optimize import minimize
#from scipy import optimize
import pickle
import numpy as np
def boostrap(seed,
             dataZ = z_data, dataX = x_data, 
             dataAiMcd = ai_mcd[tertile==mk].reshape([-1,TimePeriod]),
             dataAiOther = ai_other[tertile_other==mk].reshape([-1,TimePeriod]), 
             dataStateMcd=stateID_mcd[tertile==mk].reshape([-1,TimePeriod]),
             dataStateOther = stateID_other[tertile_other==mk].reshape([-1,TimePeriod])):
    np.random.seed(seed)
    samples = np.random.randint(0,len(dataX),len(dataX))#np.arange(len(dataX),dtype = int)#
    print(samples)
    dataZ_new = np.zeros([TimePeriod,len(dataX)],dtype = int)
    haha = TimePeriod*len(dataX)
    dataAiMcd_new= np.zeros(haha,dtype = int)
    dataAiOther_new= np.zeros(4*haha,dtype = int)
    dataStateMcd_new= np.zeros(haha,dtype = int)
    dataStateOther_new= np.zeros(4*haha,dtype = int)
    for i,ss in enumerate(samples):
        j = TimePeriod*i
        dataZ_new[:,i] = dataZ[:,ss]
        dataAiMcd_new[j:j+TimePeriod] = dataAiMcd[ss]
        dataAiOther_new[j:j+TimePeriod] = dataAiOther[ss]
        dataAiOther_new[j+haha:j+TimePeriod+haha] = dataAiOther[ss+len(dataX)]
        dataAiOther_new[j+2*haha:j+TimePeriod+2*haha] = dataAiOther[ss+2*len(dataX)]
        dataAiOther_new[j+3*haha:j+TimePeriod+3*haha] = dataAiOther[ss+3*len(dataX)]

        dataStateMcd_new[j:j+TimePeriod] = dataStateMcd[ss]
        dataStateOther_new[j:j+TimePeriod] = dataStateOther[ss]
        dataStateOther_new[j+haha:j+TimePeriod+haha] = dataStateOther[ss+len(dataX)]
        dataStateOther_new[j+2*haha:j+TimePeriod+2*haha] = dataStateOther[ss+2*len(dataX)]
        dataStateOther_new[j+3*haha:j+TimePeriod+3*haha] = dataStateOther[ss+3*len(dataX)]

    data = {'dataZ':dataZ_new,'dataAiMcd':dataAiMcd_new,'dataAiOther':dataAiOther_new,'dataStateMcd':dataStateMcd_new,'dataStateOther':dataStateOther_new}
    #dat_new = pd.DataFrame(data=data)
    #print(samples)
    #print(dat_new)
    return data

for seed in range(seed_ini,seed_ini +10):
  print(seed)
  dat_new = boostrap(seed)
  DataRes = []
  DataRes.append(seed)
  lower_bound = 1e-3
  a1 = 0.8993
  a2 = 0.008
  def func1(theta23):
    theta23 = np.array([a1,a2,#np.array([haha[2021]['month'][0],haha[2021]['month'][1],
                        theta23[0]])#,theta23[1]])
    sigma_val,_,_ = SigmaLambda(theta23,z_old = dat_new['dataZ'],x_old=x_data,T=TimePeriod)
    #sigma_val[np.isnan(sigma_val)] = 0
    res = -np.sum(np.log(sigma_val[np.where(sigma_val!=0)]))

    #print(theta23,res)
    return res /(5*MT)

  rranges = [(0.5+lower_bound,1-lower_bound)]#,(lower_bound,0.5-lower_bound)]# *2
  init_x0 = [0.9]#,0.1]#np.array([trans_gdp[0][0],trans_gdp[1,0],0.9+lower_bound])# ,0.5-lower_bound])#np.array([2.47142617e-01, 6.98155974e-02, 9.50769569e-01, 1.00000000e-06])#
  #print('start from:',init_x0)
  res1 = minimize(func1,init_x0,bounds = rranges,options={'maxiter':1000},tol=1e-5)
  print(res1,res1.fun)#*5*MT)
  print(Dynamic(np.array([a1,a2,#np.array([haha[2021]['month'][0],haha[2021]['month'][1],
                        res1.x[0]])))#,res1.x[1]])))
  # resbrute = optimize.brute(func1, rranges, finish=None,Ns=10)
  # print(resbrute[0])
  DataRes.append(res1)
  norm_thetaS = 1
  sigma1,_,(iterp1,iterf1,iterc1) = SigmaLambda(np.array([a1,a2,#np.array([haha[2021]['month'][0],haha[2021]['month'][1],
                        res1.x[0]]))
  sigmaP1 = np.kron(sigma1.reshape([NX,obser,1]),np.ones([1,1,NS]))
  iterpP1 = np.kron(iterp1.reshape([1,NX,obser,1]),np.ones([J,1,1,NS]))

  _,_,(iterp2,iterf2,iterc2) = SigmaLambda(np.array([a1,a2,#np.array([haha[2021]['month'][0],haha[2021]['month'][1],
                        res1.x[0]]),z_old = dat_new['dataZ'],x_old=x_data,T=TimePeriod)
  #x_blief = (x_blief.T).reshape(-1)
  iterp2 = (iterp2.T).reshape(-1)
  iterf2 = (iterf2.T).reshape(-1)
  iterc2 = (iterc2.T).reshape(-1)
  iterp2_other = np.hstack([iterp2 for i in range(4)])
  iterf2_other = np.hstack([iterf2 for i in range(4)])
  iterc2_other = np.hstack([iterc2 for i in range(4)])

  def func(thetaS):
    #thetaS = np.hstack([thetaS[0:7],0/norm_thetaS,thetaS[7:15],0/norm_thetaS,thetaS[15]])
    _, Pi_new,_ = ValueIteration(thetaS*norm_thetaS,sigmaP=sigmaP1,iterpP=iterpP1,iterf=iterf1,iterc=iterc1)
    #print(np.unique(a_mcd))
    res = -np.sum(np.log(  iterp2 *Pi_new[0,dat_new['dataAiMcd'],iterf2,dat_new['dataStateMcd']]+\
                        (1-iterp2)*Pi_new[0,dat_new['dataAiMcd'],iterc2,dat_new['dataStateMcd']]))
    res = res -np.sum(np.log(  iterp2_other *Pi_new[1,dat_new['dataAiOther'],iterf2_other,dat_new['dataStateOther']]+\
                            (1-iterp2_other)*Pi_new[1,dat_new['dataAiOther'],iterc2_other,dat_new['dataStateOther']]))
    print(thetaS*norm_thetaS,res)
    return res/5/MT


  thetaS = thetaS/norm_thetaS
  print('start from:',thetaS*norm_thetaS )
  res22 = minimize(func,thetaS,options={'maxiter':4000},tol=1e-3)
  print(res22,res22.fun*5*MT) # 28min

  haha = res22.x*norm_thetaS
  print('{:.4f}&{:.4f}'.format(haha[0],haha[8]))
  print('{:.4f}&{:.4f}'.format(haha[1],haha[9]))
  print('{:.4f}&{:.4f}'.format(haha[2],haha[10]))
  print('{:.4f}&{:.4f}'.format(haha[3],haha[11]))
  print('{:.4f}&{:.4f}'.format(haha[4],haha[12]))
  print('{:.4f}&{:.4f}'.format(haha[5],haha[13]))
  print('{:.4f}&{:.4f}'.format(haha[6],haha[14]))
  print('{:.4f}&{:.4f}'.format(haha[7],haha[15]))
  DataRes.append(res22)
  DataRes.append(res22)
  with open(file_name+'Full_Mk'+str(mk)+'_Seed'+str(seed)+'.txt','wb') as fp:
    pickle.dump(DataRes,fp)
    

#%%


# mat = scipy.io.loadmat(file_name + 'canadafastfood_resorted.mat')
# data = mat['data']#[np.where(mtype_int==1)]
# TimePeriod = len(np.unique(data[:,1]))
# SamplePath = int(len(data)/TimePeriod)
# print(TimePeriod,SamplePath)
# print(data[:,0].reshape([TimePeriod,SamplePath]).T)

# #% Labeling
# MT = len(data)          #% Size of the dataset (M x T = 400 x 35 = 14000).#clusterid = np.int_(data[:,0])      % Unique market ID.
# year = np.int_(data[:,1]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# N_aw  = np.int_(data[:,2]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# N_bk  = np.int_(data[:,3]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# N_hvy = np.int_(data[:,4]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# N_mcd = np.int_(data[:,5]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# N_wdy = np.int_(data[:,6]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# lagN_aw  = np.int_(data[:,7]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# lagN_bk  = np.int_(data[:,8]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# lagN_hvy = np.int_(data[:,9]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# lagN_mcd = np.int_(data[:,10]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# lagN_wdy = np.int_(data[:,11]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fwdN_aw  = np.int_(data[:,12]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fwdN_bk  = np.int_(data[:,13]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fwdN_hvy = np.int_(data[:,14]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fwdN_mcd = np.int_(data[:,15]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fwdN_wdy = np.int_(data[:,16]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# a_aw  = np.int_(data[:,17]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# a_bk  = np.int_(data[:,18]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# a_hvy = np.int_(data[:,19]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# a_mcd = np.int_(data[:,20]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# a_wdy = np.int_(data[:,21]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
# pop = data[:,23].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# val = data[:,24].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# inc = data[:,25].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fpop = data[:,27].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# finc = data[:,28].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# fval = data[:,29].reshape([TimePeriod,SamplePath]).T.reshape(-1)
# mktfe = data[:,30].reshape([TimePeriod,SamplePath]).T.reshape(-1) #% Market fixed effect estimates (from 130913_alaToivanenWaterson3_Fixedeffect3quantile.csv/dta)
# tertile = np.int_(data[:,31]).reshape([TimePeriod,SamplePath]).T.reshape(-1)       #% Market type initial guess (from 130913_alaToivanenWaterson3_Fixedeffect3quantile.csv/dta)

# def quantile(x,p):
#   x = np.sort(x)
#   n = len(x)
#   x = np.hstack([x[0], x,x[-1]])
#   i = p*n+1.5
#   iu = int(np.ceil(i))
#   il = int(np.floor(i))
#   d = i-il
#   qq = x[il-1]*(1-d)+x[iu-1]*d
#   return qq

# #% Number of own shops (in state space; capped at 3), from the perspective of each firm
# Ni_aw = (N_aw <= 3)*N_aw + (N_aw > 3)*3
# Ni_bk = (N_bk <= 3)*N_bk + (N_bk > 3)*3
# Ni_hvy = (N_hvy <= 3)*N_hvy + (N_hvy > 3)*3
# Ni_mcd = (N_mcd <= 3)*N_mcd + (N_mcd > 3)*3
# Ni_wdy = (N_wdy <= 3)*N_wdy + (N_wdy > 3)*3

# #% Number of rival shops (in state space; capped at 3), from the perspective of each firm
# Nj_aw = N_bk + N_hvy + N_mcd + N_wdy
# Nj_bk = N_aw + N_hvy + N_mcd + N_wdy
# Nj_hvy = N_aw + N_bk + N_mcd + N_wdy
# Nj_mcd = N_aw + N_bk + N_hvy + N_wdy
# Nj_wdy = N_aw + N_bk + N_hvy + N_mcd
# Nj_aw = (Nj_aw <= 3)*Nj_aw + (Nj_aw > 3)*3     
# Nj_bk = (Nj_bk <= 3)*Nj_bk + (Nj_bk > 3)*3
# Nj_hvy = (Nj_hvy <= 3)*Nj_hvy + (Nj_hvy > 3)*3
# Nj_mcd = (Nj_mcd <= 3)*Nj_mcd + (Nj_mcd > 3)*3
# Nj_wdy = (Nj_wdy <= 3)*Nj_wdy + (Nj_wdy > 3)*3

# #% Number of pop,val,inc (in state space; capped at 3)
# pop25 = quantile(pop,0.25)
# val25 = quantile(val,0.25)   
# inc25 = quantile(inc,0.25)  
# pop50 = quantile(pop,0.5) 
# val50 = quantile(val,0.5)   
# inc50 = quantile(inc,0.5)  
# pop75 = quantile(pop,0.75) 
# val75 = quantile(val,0.75)   
# inc75 = quantile(inc,0.75)   
# disc_pop = np.int_(pop > pop25) + np.int_(pop > pop50) + np.int_(pop > pop75)  
# disc_val = np.int_(val > val25) + np.int_(val > val50) + np.int_(val > val75)  
# disc_inc = np.int_(inc > inc25) + np.int_( inc > inc50) + np.int_(inc > inc75)

# #% Define state variables for each firm
# RHS_aw = np.vstack([Ni_aw, Nj_aw, disc_pop,disc_val, disc_inc]).T#;         % State for A & W.      
# RHS_bk = np.vstack([Ni_bk, Nj_bk, disc_pop,disc_val, disc_inc]).T#;         % State for Burger King.      
# RHS_hvy = np.vstack([Ni_hvy, Nj_hvy, disc_pop,disc_val, disc_inc]).T#;      % State for Harvey's.      
# RHS_mcd = np.vstack([Ni_mcd, Nj_mcd, disc_pop, disc_val,disc_inc]).T#;      % State for McDonald's.      
# RHS_wdy = np.vstack( [Ni_wdy, Nj_wdy, disc_pop,disc_val, disc_inc]).T#;


# stateID_aw = Ni_aw*4**4+Nj_aw*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
# stateID_bk = Ni_bk*4**4+Nj_bk*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
# stateID_hvy = Ni_hvy*4**4+Nj_hvy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
# stateID_mcd = Ni_mcd*4**4+Nj_mcd*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
# stateID_wdy = Ni_wdy*4**4+Nj_wdy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
# stateID_other = np.hstack([stateID_aw,stateID_bk,stateID_hvy,stateID_wdy])

# ai_aw = -1*np.int_(a_aw<0) + 0*np.int_(a_aw==0) + 1*np.int_(a_aw>0)+1
# ai_bk = -1*np.int_(a_bk<0) + 0*np.int_(a_bk==0) + 1*np.int_(a_bk>0) +1
# ai_hvy = -1*np.int_(a_hvy<0) + 0*np.int_(a_hvy==0) + 1*np.int_(a_hvy>0)+1
# ai_mcd = -1*np.int_(a_mcd<0) + 0*np.int_(a_mcd==0) + 1*np.int_(a_mcd>0)+1
# ai_wdy = -1*np.int_(a_wdy<0) + 0*np.int_(a_wdy==0) + 1*np.int_(a_wdy>0)+1
# ai_other = np.hstack([ai_aw,ai_bk,ai_hvy,ai_wdy])

# numstates = 5 # ni,nj,pop,val,inc
# states_full = np.arange(4**numstates)
# Ni = states_full//(4**(numstates-1))                       #Number of own outlets: {0,1,2,3+}.
# Nj = states_full%(4**(numstates-1))//(4**(numstates-2))    #Number of rival outlets: {0,1,2,3+}.
# dz1 = states_full%(4**(numstates-2))//(4**(numstates-3))   #Discretized population state: {0,1,2,3}.
# dz2 = states_full%(4**(numstates-3))//(4**(numstates-4))   #Discretized property value state: {0,1,2,3}.
# dz3 = states_full%(4**(numstates-4))                       #Discretized income state: {0,1,2,3}.
# print(np.sum(np.abs(Ni*4**4+Nj*4**3+dz1*4**2+dz2*4**1+dz3-states_full)))
# states = np.stack([Ni,Nj,dz1,dz2,dz3],axis=1)


# filename = file_name#'/content/drive/MyDrive/2022_SUMMER/Research_Burger/Picture/'#file_name#
# data_extra = pd.read_csv(filename+'ExtraData.csv')
# print(data_extra.columns)
# obser = 2
# # gdp = data_extra['Inflation_Rate']#['disc_inf']#['disc_une']#['Inflation_Rate']#['GDP_Grwoth_Rate']#['Unemployment_Rate']#
# inf = data_extra['InflationRateLag_wb']#['disc_inf']#['disc_une']#['disc_gdp']
# inf50 = quantile(inf,0.5)
# print(inf50)
# disc_inf = np.int_(inf <= inf50)#0)#<gdp50)#)#np.int_(gdp <0)#
# trans_inf = np.zeros([obser,obser],dtype=int)

# for i in range(34):
#   trans_inf[disc_inf[i],disc_inf[i+1]] +=1
# print(trans_inf)
# trans_inf = trans_inf/np.kron(np.ones([1,obser]),np.sum(trans_inf,axis=1).reshape([-1,1]))
# print(trans_inf)
# disc_inf = (np.repeat(disc_inf[0:TimePeriod],SamplePath).reshape([TimePeriod,SamplePath]).T).reshape(-1)


# gdp = data_extra['GDP_Grwoth_Rate']#['disc_inf']#['disc_une']#['disc_gdp']
# gdp50 = quantile(gdp,0.5)
# print(gdp50)
# disc_gdp = np.int_(gdp>0)#<gdp50)#)#np.int_(gdp <0)#
# trans_gdp = np.zeros([obser,obser],dtype=int)

# for i in range(34):
#   trans_gdp[disc_gdp[i],disc_gdp[i+1]] +=1
# print(trans_gdp)
# trans_gdp = trans_gdp/np.kron(np.ones([1,obser]),np.sum(trans_gdp,axis=1).reshape([-1,1]))
# print(trans_gdp)
# disc_gdp = (np.repeat(disc_gdp[0:TimePeriod],SamplePath).reshape([TimePeriod,SamplePath]).T).reshape(-1)

# dataFull = pd.DataFrame()
# dataFull['tertile'] = tertile
# dataFull['disc_pop'] = disc_pop
# dataFull['disc_val'] = disc_val
# dataFull['disc_inc'] = disc_inc

# dataFull['Ni_aw'] = Ni_aw
# dataFull['Ni_bk'] = Ni_bk
# dataFull['Ni_hvy'] = Ni_hvy
# dataFull['Ni_mcd'] = Ni_mcd
# dataFull['Ni_wdy'] = Ni_wdy

# dataFull['Nj_aw'] = Nj_aw
# dataFull['Nj_bk'] = Nj_bk
# dataFull['Nj_hvy'] = Nj_hvy
# dataFull['Nj_mcd'] = Nj_mcd
# dataFull['Nj_wdy'] = Nj_wdy

# dataFull['ai_aw'] = ai_aw
# dataFull['ai_bk'] = ai_bk
# dataFull['ai_hvy'] = ai_hvy
# dataFull['ai_mcd'] = ai_mcd
# dataFull['ai_wdy'] = ai_wdy

# dataFull['disc_inf'] = disc_inf
# dataFull['disc_gdp'] = disc_gdp
# dataFull.to_csv(filename+'dataFull.csv',index=False)



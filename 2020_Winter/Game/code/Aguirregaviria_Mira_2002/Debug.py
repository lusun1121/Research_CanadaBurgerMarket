# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#%%
import numpy as np
from scipy.special import logsumexp
from scipy.special import softmax
import warnings
warnings.filterwarnings("ignore")

theta1 = np.array([1.5,1,-1.9,-1.8,1]) #RS, FC1,FC2, EC

trans_true = [0.9,0.1]
obser_true = [0.2,0.3]
theta23 = trans_true + obser_true

agent = 2
za_dim =(agent+1)*(2**agent)
za_old = np.arange(za_dim,dtype = int)

TimePeriod= 150
SamplePath = 3000
SubSample = 30

# z_old = np.arange(agent+1,dtype=int) #observation
# a1_old = np.arange(2,dtype = int)
# a2_old = np.arange(2,dtype = int)

num_discrete = 50
x_old = np.linspace(0,1,num_discrete)

za_old = np.kron(np.ones(2,dtype=int),np.kron(np.ones(num_discrete,dtype = int),za_old))
x_old = np.kron(np.ones(2,dtype = int),np.kron(x_old,np.ones(za_dim,dtype=int))) #z,a2_old,a1_old,x,a_new
a_new = np.kron(np.arange(2,dtype=int),np.ones(num_discrete*za_dim,dtype=int))
#%%

def Dynamic(theta23):

  trans_val = theta23[0:2]
  obser_val = theta23[2:4]
  trans = np.zeros([2,2]) #s_t,s_{t+1}
  trans[:,0] = trans_val
  trans[:,1] = 1- trans[:,0]
  
  observe = np.zeros([2,3,3]) #s_{t+1},z_t,z_{t+1}
  for s in range(2):
    observe[s,0,1] = observe[s,1,0] = observe[s,1,2]= observe[s,2,1] = obser_val[s]
    observe[s,0,0] = 1 - observe[s,0,1] -observe[s,0,2]
    observe[s,1,1] = 1 - observe[s,1,0] -observe[s,1,2]
    observe[s,2,2] = 1 - observe[s,2,0] - observe[s,2,1]

  return trans,observe
#print(Dynamic(theta23))

def Reward(theta1,ag_i): # 3*2*2*100*2,2
  rw = np.zeros([(agent+1)*(2**agent)*num_discrete*2,2])#z*a_{t-1}*x_t*a_{t,i},a_{t,-i}
  #for s in range(2):
  if ag_i==1:
    rw[:,0] = ((za_old%(agent+1)) * theta1[0]/((2+0)**2) - theta1[2] - theta1[4] * (1-((za_old//(agent+1))//2)))*x_old
    rw[:,0] += ((za_old%(agent+1)) * theta1[1]/((2+0)**2) - theta1[2] - theta1[4] * (1-((za_old//(agent+1))//2)))*(1-x_old)

    rw[:,1] = ((za_old%(agent+1)) * theta1[0]/((2+1)**2) - theta1[2] - theta1[4] * (1-((za_old//(agent+1))//2)))*x_old
    rw[:,1] += ((za_old%(agent+1)) * theta1[1]/((2+1)**2) - theta1[2] - theta1[4] * (1-((za_old//(agent+1))//2)))*(1-x_old)

  elif ag_i==2:
    rw[:,0] = ((za_old%(agent+1)) * theta1[0]/((2+0)**2) - theta1[3] - theta1[4] * (1-((za_old//(agent+1))%2)))*x_old
    rw[:,0] += ((za_old%(agent+1)) * theta1[1]/((2+0)**2) - theta1[3] - theta1[4] * (1-((za_old//(agent+1))%2)))*(1-x_old)

    rw[:,1] = ((za_old%(agent+1)) * theta1[0]/((2+1)**2) - theta1[3] - theta1[4] * (1-((za_old//(agent+1))%2)))*x_old
    rw[:,1] += ((za_old%(agent+1)) * theta1[1]/((2+1)**2) - theta1[3] - theta1[4] * (1-((za_old//(agent+1))%2)))*(1-x_old)

  rw[:,0] = a_new * rw[:,0]
  rw[:,1] = a_new * rw[:,1]

  return rw
#print(Reward(theta1,2)[1200:2400])


def SigmaLambda(za_old,x_old,trans,observe,num_discrete=num_discrete,T = None): #(xt,z_t,a_old,a_new_i),z_{t+1}
  if T==None: #generate function for Q function iteration
      x_new = np.zeros([x_old.shape[0],3])
      sigma = np.zeros([x_old.shape[0],3])
      x_old = np.vstack([x_old,1-x_old])
      for z_new in range(3):
        x_new_tep = trans[:,0].reshape([-1,1])
        x_new_tep = x_new_tep.dot(observe[0,za_old%(agent+1),z_new].reshape([1,-1]))
        x_new_tep = np.sum(x_old * x_new_tep,axis=0)
        #print(x_new_tep.shape)
        sigma[:,z_new] = x_new_tep + np.sum(x_old*(np.reshape(trans[:,1],[-1,1]).dot(observe[1,za_old%(agent+1),z_new].reshape([1,-1]))),
                                            axis=0)
        sigma_0 = np.where(sigma[:,z_new]!=0)
        x_new[:,z_new][sigma_0] = x_new_tep[sigma_0] /sigma[:,z_new][sigma_0] 

      belief_f = np.floor(x_new*(num_discrete-1))/(num_discrete-1)
      belief_c = np.ceil(x_new*(num_discrete-1))/(num_discrete-1)
      iterpolate = np.zeros(belief_f.shape)
      iterpolate[np.where((belief_f-belief_c)!=0)] = (x_new-belief_c)[np.where((belief_f-belief_c)!=0)]/(belief_f-belief_c)[np.where((belief_f-belief_c)!=0)]

  else:      #generate whole belief for recover process
      x_new = np.zeros([T,x_old.shape[0]]) # T, N
      sigma = np.zeros([T-1,x_old.shape[0]]) # T-1 ,N
      x_new[0,:] = x_old

      for t in range(1,T):

        za = za_old[t-1,:]
        z = za_old[t,:]
        
        x_new_tep = (x_new[t-1,:] *trans[0,0] + (1-x_new[t-1,:])*trans[1,0] )* observe[0,za,z]
        sigma[t-1,:] = x_new_tep + (x_new[t-1,:] *trans[0,1] + (1-x_new[t-1,:])*trans[1,1] )* observe[1,za,z]
        x_new[t,:] = x_new_tep/sigma[t-1,:]

      belief_f = np.floor(x_new*(num_discrete-1))/(num_discrete-1)
      belief_c = np.ceil(x_new*(num_discrete-1))/(num_discrete-1)
      iterpolate = np.zeros(belief_f.shape)
      iterpolate[np.where((belief_f-belief_c)!=0)] = (x_new-belief_c)[np.where((belief_f-belief_c)!=0)]/(belief_f-belief_c)[np.where((belief_f-belief_c)!=0)]      

  return sigma,x_new,[iterpolate,np.int_(belief_f*(num_discrete-1)),np.int_(belief_c*(num_discrete-1))]
# trans,observe = Dynamic(theta23)  #s_t,s_{t+1} #s_{t+1},z_t,z_{t+1}
# print(SigmaLambda(za_old,x_old,trans,observe)[0].shape)

# trans,observe = Dynamic(theta23)  #s_t,s_{t+1} #s_{t+1},z_t,z_{t+1}
# sigma,x_new,xiterpolate = SigmaLambda(za_old,x_old,trans,observe)
#%%

myZero = 1e-5

import time
import scipy
def ValueIteration(theta1,sigma,xiterpolate, num_discrete=num_discrete,beta = 0.95,
                   gamma = 0.5772,thread_inner = 1e-5,thread_outer = 1e-5): #xt*za*a1 or xt*za*a2
  rw1,rw2 = Reward(theta1,1), Reward(theta1,2) #z*a_{t-1}*xt*a_{t,i},a_{t,-i}
  iterpolate,belief_f,belief_c = xiterpolate[0],xiterpolate[1],xiterpolate[2]#z*a_{t-1}*xt*a_{t,i},a_{t,-i}
  x_old = np.linspace(0,1,num=num_discrete)#.reshape([-1,1])
  x_old = np.kron(np.ones(2,dtype = int),np.kron(x_old,np.ones(za_dim,dtype=int))) #z,a2_old,a1_old,x,a_new  
  #x_old = np.vstack([x_old,1-x_old]).T

  pi1_new = 0.5*np.ones(num_discrete*za_dim*2) #xt*za*a1

  pi2_old = pi1_new.copy()
  #count = 0
  #pi1_old = np.zeros([2,2])
  for ol in range(1000):
    pi1_old = pi1_new.copy()
    
    #Q1_new = np.zeros(num_discrete*za_dim*2) #xt,za,a1
    h1_new = np.zeros(num_discrete*za_dim*2) #xt,za,a1
    g1_new = np.zeros(num_discrete*za_dim*2) #xt,za,a1
    
    h2_new = np.zeros(num_discrete*za_dim*2) #xt,za,a2
    g2_new = np.zeros(num_discrete*za_dim*2) #xt,za,a1


    Idmat1 = np.zeros([num_discrete*za_dim*2,num_discrete*za_dim*2])
    xx = np.arange(num_discrete*za_dim*2,dtype=int)
    for z_new in range(3):
      znew_a10= (0 *2+a_new)*3 + z_new
      it0,bf0,bc0 = iterpolate[:,z_new],belief_f[:,z_new],belief_c[:,z_new]
      #print(znew_a10,bf0)
      Idmat1[xx,bf0*za_dim+znew_a10+0*num_discrete*za_dim] +=    it0*np.kron(np.ones(2),pi2_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat1[xx,bc0*za_dim+znew_a10+0*num_discrete*za_dim] +=(1-it0)*np.kron(np.ones(2),pi2_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat1[xx,bf0*za_dim+znew_a10+1*num_discrete*za_dim] +=    it0*np.kron(np.ones(2),pi2_old[num_discrete*za_dim:])*sigma[:,z_new]
      Idmat1[xx,bc0*za_dim+znew_a10+1*num_discrete*za_dim] +=(1-it0)*np.kron(np.ones(2),pi2_old[num_discrete*za_dim:])*sigma[:,z_new]
    pi1_tep1  = np.kron(np.ones(2),pi1_old[0:num_discrete*za_dim])
    Idmat1 = np.dot(pi1_tep1.reshape([-1,1]),np.ones([1,num_discrete*za_dim*2])) * Idmat1
    
    Idmat2 = np.zeros([num_discrete*za_dim*2,num_discrete*za_dim*2])
    for z_new in range(3):
      znew_a11= (1 *2+a_new)*3 + z_new
      it1,bf1,bc1 = iterpolate[:,z_new],belief_f[:,z_new],belief_c[:,z_new]
      Idmat2[xx,bf1*za_dim+znew_a11+0*num_discrete*za_dim] +=    it1*np.kron(np.ones(2),pi2_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat2[xx,bc1*za_dim+znew_a11+0*num_discrete*za_dim] +=(1-it1)*np.kron(np.ones(2),pi2_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat2[xx,bf1*za_dim+znew_a11+1*num_discrete*za_dim] +=    it1*np.kron(np.ones(2),pi2_old[num_discrete*za_dim:])*sigma[:,z_new]
      Idmat2[xx,bc1*za_dim+znew_a11+1*num_discrete*za_dim] +=(1-it1)*np.kron(np.ones(2),pi2_old[num_discrete*za_dim:])*sigma[:,z_new]
    pi1_tep2  = np.kron(np.ones(2),pi1_old[num_discrete*za_dim:])
    Idmat2 = np.dot(pi1_tep2.reshape([-1,1]),np.ones([1,num_discrete*za_dim*2])) * Idmat2
    
    Idmat1 = beta * (Idmat1 + Idmat2)
    
    for il in range(1000):
      h2_old = h2_new.copy() #xt,za,a2
      g2_old = g2_new.copy() #xt,za,a2
      
      h2_new =   pi1_tep1 * rw2[:,0] + pi1_tep2 * rw2[:,1] + Idmat1.dot(h2_old)      
      g2_new =  Idmat1.dot(gamma-np.log(pi2_old)) + Idmat1.dot(g2_old)
      #print(it0)
      
      print(il,np.max(np.abs(np.nan_to_num((h2_new-h2_old)/h2_old))) <thread_inner ,np.max(np.abs(np.nan_to_num((g2_new-g2_old)/g2_old))))
      if np.max(np.abs(np.nan_to_num((h2_new-h2_old)/h2_old))) <thread_inner and np.max(np.abs(np.nan_to_num((g2_new-g2_old)/g2_old))) <thread_inner: 
        #print('Q2',il,np.max(np.abs(np.nan_to_num((Q2_new-Q2_old)/Q2_old))))
        break

    pi2_new = np.reshape(softmax(np.reshape(h2_new + g2_new,[2,-1]),axis=0),-1)   
    count2 = np.max(np.abs(np.nan_to_num((pi2_old-pi2_new)/pi2_old)))


    pi2_old = pi2_new.copy()
    
    Idmat1 = np.zeros([num_discrete*za_dim*2,num_discrete*za_dim*2])
    xx = np.arange(num_discrete*za_dim*2,dtype=int)
    for z_new in range(3):
      znew_a20= (a_new*2+0)*3 + z_new
      it0,bf0,bc0 = iterpolate[:,z_new],belief_f[:,z_new],belief_c[:,z_new]
      #print(znew_a10,bf0)
      Idmat1[xx,bf0*za_dim+znew_a20+0*num_discrete*za_dim] +=    it0*np.kron(np.ones(2),pi1_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat1[xx,bc0*za_dim+znew_a20+0*num_discrete*za_dim] +=(1-it0)*np.kron(np.ones(2),pi1_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat1[xx,bf0*za_dim+znew_a20+1*num_discrete*za_dim] +=    it0*np.kron(np.ones(2),pi1_old[num_discrete*za_dim:])*sigma[:,z_new]
      Idmat1[xx,bc0*za_dim+znew_a20+1*num_discrete*za_dim] +=(1-it0)*np.kron(np.ones(2),pi1_old[num_discrete*za_dim:])*sigma[:,z_new]
    pi2_tep1  = np.kron(np.ones(2),pi2_old[0:num_discrete*za_dim])
    Idmat1 = np.dot(pi2_tep1.reshape([-1,1]),np.ones([1,num_discrete*za_dim*2])) * Idmat1
    
    Idmat2 = np.zeros([num_discrete*za_dim*2,num_discrete*za_dim*2])
    for z_new in range(3):
      znew_a21= (a_new*2+1)*3 + z_new
      it1,bf1,bc1 = iterpolate[:,z_new],belief_f[:,z_new],belief_c[:,z_new]
      Idmat2[xx,bf1*za_dim+znew_a21+0*num_discrete*za_dim] +=    it1*np.kron(np.ones(2),pi1_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat2[xx,bc1*za_dim+znew_a21+0*num_discrete*za_dim] +=(1-it1)*np.kron(np.ones(2),pi1_old[0:num_discrete*za_dim])*sigma[:,z_new]
      Idmat2[xx,bf1*za_dim+znew_a21+1*num_discrete*za_dim] +=    it1*np.kron(np.ones(2),pi1_old[num_discrete*za_dim:])*sigma[:,z_new]
      Idmat2[xx,bc1*za_dim+znew_a21+1*num_discrete*za_dim] +=(1-it1)*np.kron(np.ones(2),pi1_old[num_discrete*za_dim:])*sigma[:,z_new]
    pi2_tep2  = np.kron(np.ones(2),pi2_old[num_discrete*za_dim:])
    Idmat2 = np.dot(pi2_tep2.reshape([-1,1]),np.ones([1,num_discrete*za_dim*2])) * Idmat2
    
    Idmat1 = beta * (Idmat1 + Idmat2)
    
    for il in range(1000):
      h1_old = h1_new.copy() #xt,za,a2
      g1_old = g1_new.copy() #xt,za,a2


      h1_new =   pi2_tep1 * rw1[:,0] + pi2_tep2 * rw1[:,1] + Idmat1.dot(h1_old)      
      g1_new =  Idmat1.dot(gamma-np.log(pi1_old)) + Idmat1.dot(g1_old)
      
      #print(it0)
      
      print(il,np.max(np.abs(np.nan_to_num((h1_new-h1_old)/h1_old))) <thread_inner ,np.max(np.abs(np.nan_to_num((g1_new-g1_old)/g1_old))))
      if np.max(np.abs(np.nan_to_num((h1_new-h1_old)/h1_old))) <thread_inner and np.max(np.abs(np.nan_to_num((g1_new-g1_old)/g1_old))) <thread_inner: 
        #print('Q1',il,np.max(np.abs(np.nan_to_num((Q1_new-Q1_old)/Q1_old))))
        break

    pi1_new = np.reshape(softmax(np.reshape(h1_new + g1_new,[2,-1]),axis=0),-1)   
    count1 = np.max(np.abs(np.nan_to_num((pi1_old-pi1_new)/pi1_old)))
    
    #print('pi1 pi2',ol,count1,count2)

    if count1 <thread_outer and count2<thread_outer:
      #print(ol,np.max(np.abs((pi1_old-pi1_new)/pi1_new)))
      break
  return pi1_new,pi2_new,[h1_new,g1_new],[h2_new,g2_new]

trans,observe = Dynamic(theta23)  #s_t,s_{t+1} #s_{t+1},z_t,z_{t+1}
sigma,x_new,xiterpolate = SigmaLambda(za_old,x_old,trans,observe)
pi1,pi2,Q1,Q2 = ValueIteration(theta1,sigma,xiterpolate)
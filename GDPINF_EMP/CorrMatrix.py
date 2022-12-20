import numpy as np

#file_name = 'G:/My Drive/2022_SUMMER/Research_Burger/CODE/TranslatePythonCode/'
file_name = 'D:/Google Drive/2022_SUMMER/Research_Burger/CODE/TranslatePythonCode/'
import scipy.io
mat = scipy.io.loadmat(file_name + 'canadafastfood_resorted.mat')
data = mat['data']#[np.where(mtype_int==1)]
TimePeriod = len(np.unique(data[:,1]))
SamplePath = int(len(data)/TimePeriod)
print(TimePeriod,SamplePath)
print(data[:,0].reshape([TimePeriod,SamplePath]).T)

#% Labeling
MT = len(data)          #% Size of the dataset (M x T = 400 x 35 = 14000).#clusterid = np.int_(data[:,0])      % Unique market ID.
year = np.int_(data[:,1]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
N_aw  = np.int_(data[:,2]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
N_bk  = np.int_(data[:,3]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
N_hvy = np.int_(data[:,4]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
N_mcd = np.int_(data[:,5]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
N_wdy = np.int_(data[:,6]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
lagN_aw  = np.int_(data[:,7]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
lagN_bk  = np.int_(data[:,8]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
lagN_hvy = np.int_(data[:,9]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
lagN_mcd = np.int_(data[:,10]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
lagN_wdy = np.int_(data[:,11]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
fwdN_aw  = np.int_(data[:,12]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
fwdN_bk  = np.int_(data[:,13]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
fwdN_hvy = np.int_(data[:,14]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
fwdN_mcd = np.int_(data[:,15]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
fwdN_wdy = np.int_(data[:,16]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
a_aw  = np.int_(data[:,17]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
a_bk  = np.int_(data[:,18]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
a_hvy = np.int_(data[:,19]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
a_mcd = np.int_(data[:,20]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
a_wdy = np.int_(data[:,21]).reshape([TimePeriod,SamplePath]).T.reshape(-1)
pop = data[:,23].reshape([TimePeriod,SamplePath]).T.reshape(-1)
val = data[:,24].reshape([TimePeriod,SamplePath]).T.reshape(-1)
inc = data[:,25].reshape([TimePeriod,SamplePath]).T.reshape(-1)
fpop = data[:,27].reshape([TimePeriod,SamplePath]).T.reshape(-1)
finc = data[:,28].reshape([TimePeriod,SamplePath]).T.reshape(-1)
fval = data[:,29].reshape([TimePeriod,SamplePath]).T.reshape(-1)
mktfe = data[:,30].reshape([TimePeriod,SamplePath]).T.reshape(-1) #% Market fixed effect estimates (from 130913_alaToivanenWaterson3_Fixedeffect3quantile.csv/dta)
tertile = np.int_(data[:,31]).reshape([TimePeriod,SamplePath]).T.reshape(-1)       #% Market type initial guess (from 130913_alaToivanenWaterson3_Fixedeffect3quantile.csv/dta)


def quantile(x,p):
  x = np.sort(x)
  n = len(x)
  x = np.hstack([x[0], x,x[-1]])
  i = p*n+1.5
  iu = int(np.ceil(i))
  il = int(np.floor(i))
  d = i-il
  qq = x[il-1]*(1-d)+x[iu-1]*d
  return qq

#% Number of own shops (in state space; capped at 3), from the perspective of each firm
Ni_aw = (N_aw <= 3)*N_aw + (N_aw > 3)*3
Ni_bk = (N_bk <= 3)*N_bk + (N_bk > 3)*3
Ni_hvy = (N_hvy <= 3)*N_hvy + (N_hvy > 3)*3
Ni_mcd = (N_mcd <= 3)*N_mcd + (N_mcd > 3)*3
Ni_wdy = (N_wdy <= 3)*N_wdy + (N_wdy > 3)*3

#% Number of rival shops (in state space; capped at 3), from the perspective of each firm
Nj_aw = N_bk + N_hvy + N_mcd + N_wdy
Nj_bk = N_aw + N_hvy + N_mcd + N_wdy
Nj_hvy = N_aw + N_bk + N_mcd + N_wdy
Nj_mcd = N_aw + N_bk + N_hvy + N_wdy
Nj_wdy = N_aw + N_bk + N_hvy + N_mcd
Nj_aw = (Nj_aw <= 3)*Nj_aw + (Nj_aw > 3)*3     
Nj_bk = (Nj_bk <= 3)*Nj_bk + (Nj_bk > 3)*3
Nj_hvy = (Nj_hvy <= 3)*Nj_hvy + (Nj_hvy > 3)*3
Nj_mcd = (Nj_mcd <= 3)*Nj_mcd + (Nj_mcd > 3)*3
Nj_wdy = (Nj_wdy <= 3)*Nj_wdy + (Nj_wdy > 3)*3

#% Number of pop,val,inc (in state space; capped at 3)
pop25 = quantile(pop,0.25)
val25 = quantile(val,0.25)   
inc25 = quantile(inc,0.25)  
pop50 = quantile(pop,0.5) 
val50 = quantile(val,0.5)   
inc50 = quantile(inc,0.5)  
pop75 = quantile(pop,0.75) 
val75 = quantile(val,0.75)   
inc75 = quantile(inc,0.75)   
disc_pop = np.int_(pop > pop25) + np.int_(pop > pop50) + np.int_(pop > pop75)  
disc_val = np.int_(val > val25) + np.int_(val > val50) + np.int_(val > val75)  
disc_inc = np.int_(inc > inc25) + np.int_( inc > inc50) + np.int_(inc > inc75)

#% Define state variables for each firm
RHS_aw = np.vstack([Ni_aw, Nj_aw, disc_pop,disc_val, disc_inc]).T#;         % State for A & W.      
RHS_bk = np.vstack([Ni_bk, Nj_bk, disc_pop,disc_val, disc_inc]).T#;         % State for Burger King.      
RHS_hvy = np.vstack([Ni_hvy, Nj_hvy, disc_pop,disc_val, disc_inc]).T#;      % State for Harvey's.      
RHS_mcd = np.vstack([Ni_mcd, Nj_mcd, disc_pop, disc_val,disc_inc]).T#;      % State for McDonald's.      
RHS_wdy = np.vstack( [Ni_wdy, Nj_wdy, disc_pop,disc_val, disc_inc]).T#;

stateID_aw = Ni_aw*4**4+Nj_aw*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_bk = Ni_bk*4**4+Nj_bk*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_hvy = Ni_hvy*4**4+Nj_hvy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_mcd = Ni_mcd*4**4+Nj_mcd*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_wdy = Ni_wdy*4**4+Nj_wdy*4**3+disc_pop*4**2+disc_val*4**1+disc_inc
stateID_other = np.hstack([stateID_aw,stateID_bk,stateID_hvy,stateID_wdy])

ai_aw = -1*np.int_(a_aw<0) + 0*np.int_(a_aw==0) + 1*np.int_(a_aw>0)+1
ai_bk = -1*np.int_(a_bk<0) + 0*np.int_(a_bk==0) + 1*np.int_(a_bk>0) +1
ai_hvy = -1*np.int_(a_hvy<0) + 0*np.int_(a_hvy==0) + 1*np.int_(a_hvy>0)+1
ai_mcd = -1*np.int_(a_mcd<0) + 0*np.int_(a_mcd==0) + 1*np.int_(a_mcd>0)+1
ai_wdy = -1*np.int_(a_wdy<0) + 0*np.int_(a_wdy==0) + 1*np.int_(a_wdy>0)+1
ai_other = np.hstack([ai_aw,ai_bk,ai_hvy,ai_wdy])

import pandas as pd
filename = 'D:/Google Drive/2022_SUMMER/Research_Burger/Picture/'

#filename = 'G:/My Drive/2022_SUMMER/Research_Burger/Picture/'
data_extra = pd.read_csv(filename+'ExtraData.csv')
gdp = data_extra['GDP_Grwoth_Rate']
#gdp50 = quantile(gdp,0.5)  
obser = 2
disc_gdp = np.int_(gdp > 0) #+ np.int_(gdp > gdp75)  #np.int_(gdp>gdp50)#
disc_gdp = disc_gdp
trans_gdp = np.zeros([obser,obser],dtype=int)

for i in range(34):
  trans_gdp[disc_gdp[i],disc_gdp[i+1]] +=1
print(trans_gdp)
trans_gdp = trans_gdp/np.kron(np.ones([1,obser]),np.sum(trans_gdp,axis=1).reshape([-1,1]))
print(trans_gdp)
disc_gdp = (np.repeat(disc_gdp[0:TimePeriod],SamplePath).reshape([TimePeriod,SamplePath]).T).reshape(-1)
#disc_gdp = disc_gdp[mtype_int_new==1]

import pandas as pd
data_sum = pd.DataFrame(data=[],columns=['McDonalds','A&W','Harveys','BurgerKing','Wendys','Population','Income','PropertyValue','RealGDP_wb','InflationRate','MarketType'])
data_sum['McDonalds'] = N_mcd
data_sum['A&W'] = N_aw
data_sum['Harveys'] = N_hvy
data_sum['BurgerKing'] = N_bk
data_sum['Wendys'] = N_wdy
data_sum['Population'] = pop
data_sum['Income'] = inc
data_sum['PropertyValue'] = val
data_sum['RealGDP_wb'] = np.kron(np.ones(SamplePath),data_extra['RealGDP_wb'][0:TimePeriod].to_numpy())
data_sum['InflationRate'] = np.kron(np.ones(SamplePath),data_extra['InflationRateLag_wb'][0:TimePeriod].to_numpy())
#data_sum['InflationRate'] = np.kron(np.ones(SamplePath),data_extra['InflationRate_wb'][0:TimePeriod].to_numpy())
data_sum['MarketType'] = tertile

# for i in range(1,4):
#   #np.round_(data_sum.describe().iloc[7],1)
#   print(np.round_(data_sum[data_sum['MarketType']==i].describe().iloc[7],1))
  
  
ai = np.hstack([ai_aw,ai_bk,ai_hvy,ai_mcd,ai_wdy])
ai_label=[]
for i in range(len(ai)):
  if ai[i] == 0:
    ai_label.append('exit')
  elif ai[i] ==1:
    ai_label.append('unchanged')
  else:
    ai_label.append('enter')
    
data_opr = pd.DataFrame(data=[],columns=['Action','Own','Rival','Population','Income','PropertyValue','RealGDP_wb','InflationRate','MarketType'])
data_opr['Action'] = ai_label
#data_opr['Action'] = data_opr['Action'].astype(action_type)
data_opr['Own'] = np.hstack([N_aw,N_bk,N_hvy,N_mcd,N_wdy])
data_opr['Rival'] =np.hstack([N_mcd +N_hvy +N_bk + N_wdy,
                              N_aw +N_hvy +N_mcd + N_wdy,
                              N_aw +N_mcd +N_bk + N_wdy,
                              N_aw +N_hvy +N_bk + N_wdy,
                              N_aw +N_hvy +N_bk + N_mcd])#Nj_mcd#
data_opr['Population'] = np.hstack([pop,pop,pop,pop,pop])/1000#disc_pop#
data_opr['Income'] = np.hstack([inc,inc,inc,inc,inc])/1000#disc_inc#
data_opr['PropertyValue'] = np.hstack([val,val,val,val,val])/1000#disc_val#
real_gdp = np.kron(np.ones(SamplePath),data_extra['RealGDP_wb'][0:TimePeriod].to_numpy())
data_opr['RealGDP_wb'] = np.hstack([real_gdp,real_gdp,real_gdp,real_gdp,real_gdp])/1000000#disc_gdp#
#real_inf = np.kron(np.ones(SamplePath),data_extra['InflationRate_wb'][0:TimePeriod].to_numpy())
real_inf = np.kron(np.ones(SamplePath),data_extra['InflationRateLag_wb'][0:TimePeriod].to_numpy())
data_opr['InflationRate'] = np.hstack([real_inf,real_inf,real_inf,real_inf,real_inf])#disc_infr#
data_opr['MarketType'] = np.hstack([tertile,tertile,tertile,tertile,tertile])
#%%
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#sns.set(rc={"figure.dpi":900, 'savefig.dpi':900})
data_cov = pd.DataFrame(data=[],columns=[' Real GDP',' Inflation Rate',' Income',' Population',' Property Value'])
gdp_cov = np.kron(np.ones(SamplePath),data_extra['RealGDP_wb'][0:TimePeriod].to_numpy())
#inf_cov = np.kron(np.ones(SamplePath),data_extra['InflationRate_wb'][0:TimePeriod].to_numpy())
inf_cov = np.kron(np.ones(SamplePath),data_extra['InflationRateLag_wb'][0:TimePeriod].to_numpy())
inc_cov = inc
pop_cov = pop
pro_cov = val
data_cov[' Real GDP'] = gdp_cov#(gdp_cov - min(gdp_cov))/(max(gdp_cov)-min(gdp_cov))
data_cov[' Inflation Rate'] = inf_cov#(inf_cov - min(inf_cov))/(max(inf_cov)-min(inf_cov))
data_cov[' Income'] = inc_cov#(inc_cov - min(inc_cov))/(max(inc_cov)-min(inc_cov))
data_cov[' Population'] = pop_cov#(pop_cov - min(pop_cov))/(max(pop_cov)-min(pop_cov))
data_cov[' Property Value'] = pro_cov#(pro_cov - min(pro_cov))/(max(pro_cov)-min(pro_cov))

cmap = 'coolwarm'#sns.color_palette("coolwarm")#sns.diverging_palette(220, 20, as_cmap=True)
cov_matrix = data_cov[[' Real GDP',' Income',' Population',' Property Value']].corr()#np.cov(data_cov.to_numpy().T,bias=True)
print(cov_matrix)
sns.heatmap(cov_matrix, annot=True,vmin=-1,vmax=1,fmt='.2f',cmap = cmap,linewidths=.5)
#plt.gcf().set_dpi(900)
plt.savefig('CorrGDP.png', dpi = 300,bbox_inches = 'tight')
plt.show()

cov_matrix = data_cov[[' Inflation Rate',' Income',' Population',' Property Value']].corr()#np.cov(data_cov.to_numpy().T,bias=True)
print(cov_matrix)
sns.heatmap(cov_matrix, annot=True, vmin=-1,vmax=1,fmt='.2f',cmap = cmap,linewidths=.5)
#plt.gcf().set_dpi(900)
plt.savefig('CorrINF.png', dpi = 300,bbox_inches = 'tight')
plt.show()

cov_matrix = data_cov[[' Real GDP',' Inflation Rate',' Income',' Population',' Property Value']].corr()#np.cov(data_cov.to_numpy().T,bias=True)
print(cov_matrix)
sns.heatmap(cov_matrix, annot=True, vmin=-1,vmax=1,fmt='.2f',cmap = cmap,linewidths=.5)
#plt.gcf().set_dpi(900)
plt.savefig('CorrelationMatrix.png', dpi = 300,bbox_inches = 'tight')
plt.show()


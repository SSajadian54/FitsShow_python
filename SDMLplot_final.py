import time
import numpy as np 
import pylab as py 
import matplotlib
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns
from matplotlib import rcParams
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"

from sklearn import tree
from sklearn.datasets import load_boston
from sklearn.tree import DecisionTreeRegressor
from sklearn import metrics
from sklearn import model_selection
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.model_selection import cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.tree import export_text
cmap=plt.get_cmap('viridis')
###=============================================================================
ng=int(15)
nam=["$counter$",r"$type$","$Ml$", r"$Dl$",r"$Ds$",r"$tE$", "Vt", "MaplI", "u0", "blendI", "magbI" ,r"$rho$", r"$RE(AU)$","Thre","dchi"]
f0=open("./files/distribution/dist_SDML.txt","r")
nm= sum(1 for line in f0) 
par=np.zeros((nm,ng)) 
par= np.loadtxt("./files/distribution/dist_SDML.txt")

nma=0.0
nbd=0.0;
npl=0.0
tE1=[]
tE2=[]
tE3=[]

for i in range(nm):
    count, ltyp, Ml, Dl, Ds, tE=  par[i,0], par[i,1], par[i,2], par[i,3], par[i,4], par[i,5]
    Vt, u0, blendI, magbI, ros=   par[i,6], par[i,7], par[i,8], par[i,9],par[i,10]
    RE, thre, dchi, ndet =        par[i,11],par[i,12], par[i,13],par[i,14]


    if(ltyp==10):     
        npl+=1.0
        tE1.append(tE)
    elif(ltyp==9):    
        nbd+=1.0
        tE2.append(tE)
    else:             
        nma+=1.0
        tE3.append(tE)
    
    if(tE>12.5 or ltyp<0.0 or ltyp>10.1 or Ml<0.0 or thre<0.9 or dchi<800.0 or ros<0.0 or u0>1.0001 or u0<0.0 or Ml>1.0 or ndet<3): 
        print ("There is an error !!!",    par[i,:],  count  )
        print (tE,   ltyp, Ml, thre,   dchi,    ros,   u0,    Ml, ndet)
        input("Enter a number   ")
    count=int(count) 
    '''
    if(count%10000==0):     
        f1=open("./files/distribution/l_{0:d}.dat".format(int(count)),"r")
        f2=open("./files/distribution/d_{0:d}.dat".format(int(count)),"r")
        nl=sum(1 for line in f1)
        nd=sum(1 for line in f2)
        print(nl, nd)
    #except: 
    #    nl=nd=nv=0    
    #if(nl>0 and nd>0):
        light=np.zeros((nl,2)); 
        light=np.loadtxt("./files/distribution/l_{0:d}.dat".format(count)) 
        datob=np.zeros((nd,3)); 
        datob=np.loadtxt("./files/distribution/d_{0:d}.dat".format(count)) 
        #####################################################################
        plt.clf()
        plt.cla()
        fig= plt.figure(figsize=(8,6)) 
        plt.plot(light[:,0], light[:,1], "r-",  lw= 1.5)
        plt.errorbar(datob[:,0], datob[:,1], yerr=datob[:,2], fmt='mo',capsize=0, ms=1.9, lw=1.0) 
        plt.xlabel("Time(days)")
        plt.ylabel("magnification")
        plt.title(str(round(dchi,2))  )
        fig=plt.gcf()
        fig.savefig("./histo/li_{0:d}.jpg".format(count),dpi=200)
        print("One linght is plotted!!***************",  count,  nl, nd, ndet)
    '''    
#########################################################################
print (nma*100.0/nm,   nbd*100.0/nm,    npl*100.0/nm  ,  nm,    nma+nbd+npl ) 


plt.clf()
plt.cla()
fig= plt.figure(figsize=(8,6)) 
ax1= plt.gca()
plt.hist(tE1,27, histtype='bar',ec='darkgreen',facecolor='green',alpha=1.0, rwidth=1.5, label=r"$\rm{Planets}$")
plt.hist(tE2,27, histtype='bar',ec='darkred',facecolor='red',alpha=0.5, rwidth=1.5, label=r"$\rm{Brown}~\rm{dwarfs}$")
plt.hist(tE3,27, histtype='bar',ec='darkblue',facecolor='blue',alpha=0.4, rwidth=1.5, label=r"$\rm{Main}-\rm{sequence}~\rm{stars}$")
#plt.hist(tE1+tE,27, histtype='bar',ec='darkblue',facecolor='blue',alpha=0.4, rwidth=1.5, label=r"$\rm{Main}-\rm{sequence}~\rm{stars}$")
ax1.set_ylabel(r"$\rm{Normalized}~\rm{distribution}$",fontsize=19,labelpad=0.1)
#ax1.set_xlabel(str(nam2[i]),fontsize=19, labelpad=0.1)
y_vals = ax1.get_yticks()
ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/len(tE1))) for x in y_vals]) 
y_vals = ax1.get_yticks()
plt.ylim([np.min(y_vals), np.max(y_vals)])
plt.xlim([0.2,12.0])
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.xlabel(r"$t_{\rm E}(\rm{days})$", fontsize=18)

#plt.axvline(x=np.mean(array0) , color='darkgreen', linestyle='--', lw=1.5)
#plt.axvline(x=np.mean(array1) , color='k',         linestyle='--', lw=1.5)
#plt.grid("True")
#plt.grid(linestyle='dashed')
plt.legend()
plt.legend(loc='best',fancybox=True, shadow=True)
plt.legend(prop={"size":18})
fig3= plt.gcf()
fig3.savefig("./tEs.jpg", dpi=200)


   
 

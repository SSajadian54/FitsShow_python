import numpy as np 
import matplotlib.pyplot as plt
from astropy.io import fits 
from matplotlib import gridspec
from astropy.utils.data import get_pkg_data_filename
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cm as cm
from matplotlib.colors import Normalize






###########################################################
titl=r"$\rm{MOA}-2011-\rm{BLG}-160, 2011$"

ycs,  xcs=  225.352,   512.0- 246.818  ##(reverse of nst.1 file) +the second one-512      250.696   

yps,  xps= -1.0*float(27.0+54.0/60.0+50.43/3600.0), float(17.0+58.0/60.0+38.71/3600.0)*15.0 ## DEC,  RA [degree]


pixel=float(45.0/512.0)/60.0  ## arcmin

###########################################################

def tick_fun1(xx):
    return (-(xx-xcs)*pixel+0.0)#xps) 

def tick_fun2(xx):
    return ((xx-ycs)*pixel+0.0)##yps)    
###########################################################
xx=np.arange(20 , 512.0,150.0)
rad=60

x1, x2=  int(xcs -rad/2) ,  int(xcs +rad/2)

y1, y2=  int(ycs -rad/2) ,  int(ycs +rad/2)


dat1=np.zeros((512,512))
sub= np.zeros((512,512))
sub2= np.zeros((512,512))
dat2=np.zeros((rad,rad))
dat3=np.zeros((rad,rad))
dat4=np.zeros((rad,rad))


refa=get_pkg_data_filename('./ref12a.fits') 
refb=get_pkg_data_filename('./ref12a.fits.sub.1.fits')
refc=get_pkg_data_filename('./ref12a.fits.sub.2.fits')
dat =fits.getdata(refa) 
sub =fits.getdata(refb)
sub2=fits.getdata(refc)

dat1=np.rot90(np.array(dat, dtype=float), 1)
sub= np.rot90(np.array(sub, dtype=float), 1)
sub2=np.rot90(np.array(sub2,dtype=float), 1)


for i in range(512): 
    for j in range(512): 
        if(dat1[i,j]<0.0  or dat1[i,j]==0.0): 
            dat1[i,j]=50.0
        if( sub[i,j]<0.0 or  sub[i,j]==0.0): 
            sub[i,j]=50.0  
        if( sub2[i,j]<700.0 or  sub2[i,j]==0.0): 
            sub2[i,j]=850.0    
            

dat2=dat1[x1:x2, y1:y2]
dat3=np.array(sub[x1:x2, y1:y2], dtype=float)
dat4=np.array(sub2[x1:x2,y1:y2],dtype=float)



max1=np.max(dat2);  max2=np.max(dat4);  maxm=np.max(np.array(max1, max2))
min1=np.min(dat2);  min2=np.min(dat4);  minm=np.min(np.array(min1, min2))
step1=(np.max(np.log10(dat1))- np.min(np.log10(dat1)) )/5.0
step2=abs(np.log10(maxm)-np.log10(minm))/4.0
vv=[];  uu=[]
for i in range(6): 
    vv.append(round(np.min(np.log10(dat1))+i*step1,1))
    if(i<4): uu.append(round(np.log10(minm) + i*step2 , 1 ))
    

##############################################################################
plt.close('all')
plt.clf()
fig= plt.figure(figsize=(12,9))
gs = gridspec.GridSpec(ncols=4,nrows=3)
ax0=plt.subplot(gs[:,:3])
pos=ax0.imshow(np.log10(dat1),extent=(0.0, 512, 0.0, 512.0),cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.gca().set_aspect('equal', adjustable='box')

ax0.annotate('source star',
            xy=(ycs, xcs),size=21, xycoords='data',
            xytext=(-15, -25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.4),
            horizontalalignment='right', verticalalignment='top')



ax0.set_ylim([0.0, 512.0])
ax0.set_yticks(xx)
ax0.set_yticklabels(["%.1f" % i for i in tick_fun1(xx)])
ax0.tick_params(axis="y", direction="out", length=4, width=2, color="k", labelsize=19)
plt.yticks(fontsize=23, rotation=60)
ax0.set_ylabel(r"$\rm{RA(arcm)}+$"+str(round(xps,1))+r"$^{\circ}$",fontsize=23, labelpad=0.1)

ax0.set_xlim(0.0,512.0)
ax0.set_xticks(xx)
ax0.set_xticklabels(["%.1f" % i for i in tick_fun2(xx)])
ax0.tick_params(axis="x", direction="out", length=4, width=2, color="k", labelsize=19)
plt.xticks(fontsize=23, rotation=0)
plt.subplots_adjust(wspace=0.0, hspace=0.0)
ax0.set_xlabel(r"$\rm{DEC(arcm)}+$"+str(round(yps,1))+r"$^{\circ}$",fontsize=23, labelpad=0.1)


ax0.set_title(titl,fontsize=23)
plt.subplots_adjust(hspace=.0)

axins1 = inset_axes(ax0,width="90%",height="3%", loc="upper center")
axins1.xaxis.set_ticks_position("bottom")
cb=plt.colorbar(pos,cax=axins1,orientation='horizontal', ticks=vv)
cb.set_label(label=r"$\log_{10}[\rm{count}]$", size=21, weight='bold')
cb.ax.tick_params(labelsize=21)


####################################################################

cmap='viridis'
normalizer=Normalize(uu[0],uu[len(uu)-1])
im=cm.ScalarMappable(norm=normalizer, cmap=cmap)
for i in range(3): 
    if(i==0): 
        ax1= plt.subplot(gs[2,3])
        ax1.imshow(np.log10(dat4),extent=(0, rad, 0, rad), cmap=cmap,norm=normalizer,aspect=1, interpolation='nearest', origin='lower')
        ax1.annotate(' ',xy=(rad/2,rad/2), xycoords='data',xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.05),horizontalalignment='right', verticalalignment='bottom')    
        
    if(i==1): 
        ax1= plt.subplot(gs[1,3])
        ax1.imshow(np.log10(dat3),extent=(0, rad, 0, rad), cmap=cmap,norm=normalizer, aspect=1, interpolation='nearest', origin='lower')
        ax1.annotate(' ',xy=(rad/2,rad/2), xycoords='data',xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.05),horizontalalignment='right', verticalalignment='bottom')
        
        #ax1.annotate(' ',xy=(rad/2-10,rad/2+8), xycoords='data',xytext=(-15,-25), textcoords='offset points',
        #    arrowprops=dict(facecolor='blue', shrink=0.005, width=0.05),horizontalalignment='right', verticalalignment='top')    
            
    if(i==2): 
        ax1= plt.subplot(gs[0,3])
        ax1.imshow(np.log10(dat2),extent=(0, rad, 0, rad), cmap=cmap,norm=normalizer, aspect=1, interpolation='nearest', origin='lower')
        ax1.annotate(' ',xy=(rad/2,rad/2), xycoords='data',xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.05),horizontalalignment='right', verticalalignment='bottom')
            
            
    ax1.xaxis.set_ticks_position('none') 
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    
axins1 = inset_axes(ax1,width="90%", height="3%",loc="upper center")
axins1.xaxis.set_ticks_position("bottom")
cb2=plt.colorbar(im, cax=axins1,orientation='horizontal',ticks=uu)
cb2.ax.tick_params(labelsize=23)
plt.subplots_adjust(wspace=0.0, hspace=0.0)
fig4=plt.gcf()
fig4.savefig("0160_2011.jpg", dpi=200)
####################################################################



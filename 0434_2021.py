import numpy as np 
import matplotlib.pyplot as plt
import pylab as py 
import pyfits
from astropy.io import fits 
import os, re
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from scipy.integrate import trapz
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import matplotlib as mpl
from matplotlib import gridspec
import matplotlib.dates as mdates
from matplotlib.transforms import Transform
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import(AutoLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, LogFormatterSciNotation as LogFormatter
import matplotlib.colors as colors
import matplotlib.ticker as mticker
import sys
#sys.path.append("/home/sajadian/Documents/iraf-main/")

###########################################################
def tick_fun1(xx):
    return (-(xx-xcs)*pixel+xps) 
def tick_fun2(xx):
    return ((xx-ycs)*pixel+yps)    
###########################################################
titl=r"$\rm{OGLE}-2011-\rm{BLG}-0434, 2021$"
year=int(2021)
ycs,  xcs= 225.826 , 305.097 
yps,  xps= -1.0*float(27.0+33.0/60.0+12.0/3600.0), float(18.0+3.0/60.0+26.96/3600.0)*15.0 ## DEC,  RA both in degree
pixel=float(45.0/512.0)
xx=np.arange(20 , 512.0,150.0)
rad=60
x1, x2=  int(xcs -rad/2) ,  int(xcs +rad/2)
y1, y2=  int(ycs -rad/2) ,  int(ycs +rad/2)
print "x1, x2, y1, t2:  ",   x1, x2, y1, y2

dat1=np.zeros((512,512))
dat2=np.zeros((rad,rad))
dat3=np.zeros((rad,rad))


refa=get_pkg_data_filename('./OGLE0434_21/ref21a.fits') 
refb=get_pkg_data_filename('./OGLE0434_21/ref21a.fits.sub.1.fits')
dat =fits.getdata(refa) 
sub =fits.getdata(refb)

dat1=np.array(dat, dtype=float)
dat2=dat1[x1:x2, y1:y2]
dat3=np.array(sub[x1:x2, y1:y2], dtype=float)

print dat1[int(xcs), int(ycs)],  
#input("Enter a number ")

step1=(np.max(np.log10(dat1))- np.min(np.log10(dat1)) )/5.0
step2=(np.max(np.log10(dat2))- np.min(np.log10(dat2)) )/5.0
vv=[];  uu=[]
for i in range(6): 
    vv.append(round(np.min(np.log10(dat1))+i*step1,1))
    uu.append(round(np.min(np.log10(dat2))+i*step2,1))


##############################################################################
#plt.rcParams["figure.autolayout"] = True

plt.close('all')
plt.clf()
fig= plt.figure(figsize=(9,6))
gs = gridspec.GridSpec(ncols=3,nrows=2)
ax0=plt.subplot(gs[:,:2])
pos=ax0.imshow(np.log10(dat1),extent=(0.0, 512, 0.0, 512.0),cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
plt.gca().set_aspect('equal', adjustable='box')

ax0.annotate('source star',
            xy=(ycs, xcs),size=16, xycoords='data',
            xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.4),
            horizontalalignment='right', verticalalignment='bottom')



ax0.set_ylim([0.0, 512.0])
ax0.set_yticks(xx)
ax0.set_yticklabels(["%.1f" % i for i in tick_fun1(xx)])
ax0.tick_params(axis="y", direction="out", length=4, width=2, color="k", labelsize=15)
plt.yticks(fontsize=18, rotation=20)

ax0.set_xlim(0.0,512.0)
ax0.set_xticks(xx)
ax0.set_xticklabels(["%.1f" % i for i in tick_fun2(xx)])
ax0.tick_params(axis="x", direction="out", length=4, width=2, color="k", labelsize=15)
plt.xticks(fontsize=18, rotation=0)
plt.subplots_adjust(wspace=0.0, hspace=0.0)

ax0.set_ylabel(r"$\rm{RA(arcs)}$",fontsize=18, labelpad=0.1)
ax0.set_xlabel(r"$\rm{DEC(arcs)}$",fontsize=18, labelpad=0.1)
ax0.set_title(titl,fontsize=18)
plt.subplots_adjust(hspace=.0)

axins1 = inset_axes(ax0,width="95%",height="3%", loc="upper center")
axins1.xaxis.set_ticks_position("bottom")
cb=plt.colorbar(pos,cax=axins1,orientation='horizontal', ticks=vv)
cb.set_label(label=r"$\log_{10}[\rm{count}]$", size=18, weight='bold')
cb.ax.tick_params(labelsize=17)


####################################################################

axs=[[], []]
for i in range(2): 
    if(i==0): 
        ax1= plt.subplot(gs[1,2])
        pos=ax1.imshow(np.log10(dat3),extent=(0, rad, 0, rad), cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
        #plt.plot(0.0,0.0,"*",c="k",markersize=10)
        
        ax1.annotate(' ',
            xy=(rad/2,rad/2), xycoords='data',
            xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.2),
            horizontalalignment='right', verticalalignment='bottom')

        axs[i].append(ax1)
    if(i==1): 
        ax1= plt.subplot(gs[0,2])
        ax1.imshow(np.log10(dat2),extent=(0, rad, 0, rad), cmap='viridis', aspect=1, interpolation='nearest', origin='lower')
        ax1.annotate(' ',
            xy=(rad/2,rad/2), xycoords='data',
            xytext=(-15, 25), textcoords='offset points',
            arrowprops=dict(facecolor='black', shrink=0.005, width=0.2),
            horizontalalignment='right', verticalalignment='bottom')

        axs[i].append(ax1)
    ax1.xaxis.set_ticks_position('none') 
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    
axins1 = inset_axes(ax1,width="95%", height="3%",loc="upper center")
axins1.xaxis.set_ticks_position("bottom")
cb2=plt.colorbar(pos, cax=axins1,orientation='horizontal',ticks=uu)
cb2.ax.tick_params(labelsize=17)
plt.subplots_adjust(hspace=.0)
fig4=plt.gcf()
fig4.savefig("E04342021.jpg", dpi=200)






import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
# Equinox https://data.giss.nasa.gov/cgi-bin/ar5/srevents.cgi 
# time calcu:  https://www.timeanddate.com/date/durationresult.html?d1=20&m1=3&y1=2011&d2=17&m2=5&y2=2021&h1=23&i1=31&s1=0&h2=6&i2=49&s2=47
####################################################################
##  coordinate of 0433 from Gaia (RA, DEC) 270.2963333    -27.8877500   
## Galactic coordinate (l,b)   2.640526499616078,-2.4173576746201113

ne=18; 
Dis= np.array([0.25, 0.75, 1.25, 1.75 ,2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75])#18
Ext= np.array([0.384191, 0.50733, 0.630468, 0.753606, 0.876745 , 0.999883,  1.12302, 1.24616, 1.33777,  1.39934, 1.43136, 1.45057, 1.45944, 1.46387, 1.46781, 1.46795, 1.4683, 1.4686])##18  updated  25 July 
def extinction(dl):
    ext=-1.0;   
    if(dl<Dis[0]):  
        ext= float(Ext[0] + float(Ext[0]/Dis[0])*(dl-Dis[0]))
    else:     
        for k in range(ne-1):  
            if(float((dl-Dis[k])*(dl-Dis[k+1]))<0.0 or dl==Dis[k]):  
                shib= float(Ext[k+1]-Ext[k])/(Dis[k+1]-Dis[k]) 
                ext=  float(Ext[k]+ shib*(dl-Dis[k]))
                break;   
    if(ext<0.0 or dl>Dis[ne-1]): 
        print "Error,  ext: ", ext,  "Dl:  ", dl
        input("Enter a number  ")              
    return(ext); 
###################################################################

nb=38
masm= np.zeros((nb, 2))
masm= np.loadtxt("./massm_kroupa.dat") #np.loadtxt("./FstarT_Bessel_c.dat")
def Absmagmass(mass): 
    MI=-100.0
    if(mass<masm[0,0]):      
        shib=float(masm[1,1]- masm[0,1])/(masm[1,0]- masm[0,0]) 
        MI=float(masm[0,1]+ shib*(mass-masm[0,0]) )
        
    elif(mass>masm[nb-1,0] or mass==masm[nb-1,0]):  
        shib=float(masm[nb-1,1]- masm[nb-2,1])/(masm[nb-1,0]- masm[nb-2,0]) 
        MI=float(masm[nb-1,1]+ shib*(mass-masm[nb-1,0]) )
    else: 
        for k in range(nb-1):  
            if(float((mass-masm[k,0])*(mass-masm[k+1,0]))<0.0 or  mass==masm[k,0] ):  
                shib=float(masm[k+1,1]- masm[k,1])/(masm[k+1,0]- masm[k,0]) 
                MI=float(masm[k,1]+ shib*(mass-masm[k,0]) )
                break;         
    if(MI==-100.0 or MI<-5.0 or MI>15.0):  
        print "Error:  Mass:  ", mass,  "  MI:  ",  MI
        input("Enter a number ")
    return(MI);   
#######################################3    
## 2022
pixel=float(45.0/512.0)
sou= [212.337, 304.251]#   
comp=[206.263,  313.496] 
lens=[ 221.263, 296.496]#  19.94      
mcomp= 17.408;  msou= 17.038;  mlens= 19.94  
print "2022:  Lens_source_distance:  ",  np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
ds1=np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
print "Lens_comp_distance:  ",    np.sqrt((lens[0]-comp[0])**2.0  +( lens[1]-comp[1])**2.0 )
print "source_comp_distance:  ",  np.sqrt((sou[0]-comp[0])**2.0  +( sou[1]-comp[1])**2.0 )
print  "Delta_m1:  ",     abs(msou-mlens),  "Delta_m2:  ",  abs(mlens-mcomp)
#input ("Enter a number ")    
print  "*************************************"

##2021
sou= [266.604,247.675] 
comp=[260.560,257.055]  
lens=[276.504,239.003]  
mcomp= 15.068;  msou= 14.696;  mlens= 17.680    
print "2021:  Lens_source_distance:  ",  np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
ds2=np.sqrt((lens[0]-sou[0])**2.0  +(lens[1]-sou[1])**2.0 )
print "Lens_comp_distance:  ",    np.sqrt((lens[0]-comp[0])**2.0  +( lens[1]-comp[1])**2.0 )
print "source_comp_distance:  ",  np.sqrt((sou[0]-comp[0])**2.0  +( sou[1]-comp[1])**2.0 )
print  "Delta_m1:  ",     abs(msou-mlens),  "Delta_m2:  ",  abs(mlens-mcomp)

print "delta_dec(arcs):  ", abs(lens[0]-sou[0])*pixel,  "delta_ra(arcs):  ",  abs(lens[1]-sou[1])*pixel

input("Enter a number ")

#print "<v>: ",  (ds1-ds2)*pixel*1000.0
#input ("Enter a number ")    

    
###################Constant  ######################################
year=float(365.2422)
g=0.086
dx= 0.289
pixel=float(45.0/512.0)
parcs= 30856775814913673.0
G=6.67384*pow(10.,-11.0)#;// in [m^3/s^2*kg].
velocity= 3.0*pow(10.0,8.0)#;//velosity of light
Msun=1.98892*pow(10.,30)#; //in [kg].
au= float(1.495978707*pow(10.0,11.0) )## meter 
cons= np.sqrt(4.0*G*Msun/(velocity * velocity * parcs * 1000.0) )
ti=np.zeros((2));


####################################OGLE-2011-BLG-0433 #####################################
fwhm= 5.9   ### updated 25 July
Ftot=279851.875## 7761073.00  ????
dix= float(276.504-266.604)##pixel  North-Sourth (declination)  ----  25 July updated 
diy= float(-1.0)*(239.003-247.675)##pixel East-West (right ascention)
ti[0]=float(75.0*24.0*3600.0+ 18.0*3600.0 + 6.0*60.0 + 12.0) ; ## t0 time with respect to Exoinox 2011
ti[1]=float(3710.0*24.0*3600.0+ 7.0*3600.0 + 18.0*60.0 +47.0);## 2021 image time  with .. ... 
mus_x, mus_y= -8.189, -1.085 ## Gaia data
sig_musx, sig_musy= 0.302,  0.422
parallax= 3.755718544370173 -1.4## mas
sig_par= 0.3779195  #mas
ms= 15.798 ## I-band apparement magnitude of source 
sigms= 0.000     
dif= 2.984   ### updated 25 July 
sigdif= 0.037 
tE0= 5.953; 
sig_tE=0.007
#####################################################################3


Ds0=1.0/parallax 
sig_Ds= abs(sig_par/parallax/parallax)## kpc
ml= ms+dif;
sigml= np.sqrt( sigdif*sigdif + sigms*sigms)
distp=np.sqrt(dix*dix + diy*diy )## 13.156 pixel
sigma= np.sqrt(fwhm*fwhm/(g*Ftot * 2.0*np.log(2.0)) + dx*dx *np.log(2.0)/(2.0*np.pi* fwhm*fwhm))##pixel
sig_dp= sigma*np.sqrt(2.0)#pixel
dist=distp*pixel*1000.0 #mili arcs
sig_d= sig_dp*pixel*1000.0 #mili arcs

print( "erro:dis:  ",   sig_d)

tim = float((ti[1]-ti[0])/3600.0/24.0/year)## total time in year
sig_t=float(1.0/year) ## error in time in day
mu= float(dist/tim)### mili arcs/year
sig_mu= mu* np.sqrt((sig_dp/distp)**2.0 + (sig_t/tim)**2.0)
print "Source_lens relative velocity(km/s):  ",  float(mu*0.1*au/(1000.0*year*24.0*3600.0) )
##########################################################################

ns=1300 ##1300
Ml=np.zeros((ns));  Dl= np.zeros((ns))
for i in range(ns): 
    Ml[i]=pow(10.0,-1.11 + i*1.11/ns ) #
    Dl[i]=float(0.0005 +(Ds0-0.001)*i/ns) #
###########################################################################
dl1=[];  dl2=[];  dl3=[]
ml1=[];  ml2=[];  ml3=[]
dl4=[] ; ml4=[]

for k in range(3):  
    Ds=Ds0 + float(k-1.0)*sig_Ds
    for i in range(ns):##Dl
        dl=Dl[i]
        if(dl<Ds): 
            for j in range(ns):   
                mass=Ml[j]
                xrel=float(dl/Ds)
                tetE=cons*np.sqrt(mass*(1.0-xrel)/xrel/Ds) * 180.0/np.pi * 3600.0 *1000.0 ## mili arcs
                tE= tetE*year/mu #days
                ext= extinction(dl);
                MI= Absmagmass(mass);    
                mapl=ext+5.0*np.log10(dl*100.0)+ MI   
                if(abs(tE-tE0)<sig_tE and k==1): 
                    dl1.append(dl);    ml1.append(mass); 
                if(abs(tE-tE0)<sig_tE and k==0): 
                    dl2.append(dl);    ml2.append(mass); 
                if(abs(tE-tE0)<sig_tE and k==2): 
                    dl3.append(dl);    ml3.append(mass);          
                if(abs(mapl-ml)< sigml):  
                    dl4.append(dl);    ml4.append(mass)
############################################################################
if(len(dl1)<len(dl4)): num=len(dl1);
else :                 num=len(dl4)        
tEE=np.zeros((num,4))            
tEE[:, 0]=dl1[:num] 
tEE[:, 1]=ml1[:num]    
tEE[:, 2]=dl4[:num]           
tEE[:, 3]=ml4[:num]    
fil2=open('./file0433.dat',"w"); 
np.savetxt(fil2,tEE.reshape((-1,4)) ,fmt="%.6f     %.6f    %.6f      %.6f")  
fil2.close()  

                        
############################################################################       
dd=0.0;  mm=0.0; count=0.00000023 
for i in range(len(dl1)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl1[i]-dl4[j])**2.0 + (ml1[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl1[i];   mm+=ml1[i];  count+=1.0 
DL= float(dd/count)  
ML= float(mm/count)                    
print "<DL>:   ",  DL,  "  <ML>:    ",  ML          
############################################# 
dd=0.0;  mm=0.0; count=0.0000000451
for i in range(len(dl2)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl2[i]-dl4[j])**2.0 + (ml2[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl2[i];   mm+=ml2[i];  count+=1.0 
DL1= float(dd/count)  
ML1= float(mm/count)                    
print "<DL>_Down:   ",  DL1,  "  <ML>_Down:    ",  ML1         
############################################# 
dd=0.0;  mm=0.0; count=0.000004134
for i in range(len(dl3)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl3[i]-dl4[j])**2.0 + (ml3[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl3[i];   mm+=ml3[i];  count+=1.0 
DL2= float(dd/count)  
ML2= float(mm/count)                    
print "<DL>_Upper:   ",  DL2,  "  <ML>_Upper:    ",  ML2      
           

############################################################################       
xrel=float(DL/Ds0)+ 0.00000000065 
tetE=cons*np.sqrt(ML*(1.0-xrel)/xrel/Ds0) * 180.0/np.pi * 3600.0 *1000.0 ## mili arcs
#sigDL=abs( DL * np.log(10.0)*0.4*sigml); 
#sigMass=ML* np.sqrt((sig_tE*2.0/tE0)**2.0 + (sig_mu*2.0/mu)**2.0 + (sigDL/DL/(1.0-xrel))**2.0 );
sigtetE=tetE*( sig_tE/tE0 + sig_mu/mu)  
mus= np.sqrt(mus_x*mus_x + mus_y*mus_y)
sig_mus= np.sqrt(sig_musx**2.0 + sig_musy**2.0)
param= np.array([tE0,sig_tE , mu, sig_mu , ml,sigml, DL, DL1, DL2 , ML, ML1, ML2,  tetE, sigtetE, mus, sig_mus ,Ds0, sig_Ds])
filf=open('./Tab_paper.dat',"a"); 
np.savetxt(filf,param.reshape((-1,18)), fmt="$OGLE0433$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f \pm %.2f$ & $%.2f\pm%.2f \pm %.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$  & $%.2f\pm%.2f$")
filf.close()       
                  




print "******************  OGLE-2011-BLG-0433 **********************"
print "Ds:  ", Ds0 ,  " sig_Ds:  ", sig_Ds
print "sigDL:  ", DL1,  DL2 ,    " sihmass:  ", DL1, DL2
print "mu(marcs/yrs):  ", mu, "sig_mu: ",  sig_mu 
print  "M_app(lens)[mag]:    ",  ml , "Error_Map_lens:  ",  sigml 
print "tE(days),  ", tE0, "   sig_tE: ", sig_tE
print "Dl(kpc):  ",  DL,  "  Mass(Solar mass):  ",  ML           
print "tetE (miliarcs):  ",  tetE
print "**************************************************************"


###############################################################
plt.clf()   
fig, ax = plt.subplots(figsize=(8,6))                 
plt.plot(dl1, ml1, "ro", ms=0.3)
plt.plot(dl2, ml2, "mo", ms=0.3)
plt.plot(dl3, ml3, "mo", ms=0.3)
plt.plot(dl4, ml4, "bo", ms=0.3)
ax.hlines(y=ML, xmin=0.0, xmax=DL, linestyles="--", colors="k", linewidth=1.9)
ax.vlines(x=DL, ymin=0.0, ymax=ML, linestyles="--", colors="k", linewidth=1.9)
plt.xlabel(r"$\rm{D_{l}}(kpc)$", fontsize=18)
plt.ylabel(r"$\rm{Lens}~\rm{mass}(M_{\odot})$", fontsize=18)
plt.xlim([np.min(Dl),Ds0])
plt.ylim([np.min(Ml),0.2])##np.max(Ml)])
plt.title(r"OGLE-2011-BLG-0433", fontsize=18)
plt.text(0.01,0.18,r"$t_{\rm E}(\rm{days})=$"+str(round(tE0,2))+ r"$\pm$" +str(round(sig_tE,2)), fontsize=19, color="r")
plt.text(0.01,0.16,r"$m_{\rm I,~\rm{lens}}(\rm{mag})=$"+str(round(ml,2))+r"$\pm$"+str(round(sigml,2)),fontsize=19, color="b")
plt.text(0.01,0.14,r"$D_{\rm s}(\rm{kpc})=$"+str(round(Ds0,2)),fontsize=19, color="k")
fig=plt.gcf()
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig.savefig("0433eventd.jpg",dpi=200)
py.clf()
###############################################################    
        

























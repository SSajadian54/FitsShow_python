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
####################################################################
## (l,b)= (0.14678734526340925,-1.4987587693639373)
##   RA, DEC=  267.9639583    -29.5827778 

ne=18; 
Dis= np.array([0.25, 0.75, 1.25, 1.75 ,2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75, 6.25, 6.75, 7.25, 7.75, 8.25, 8.75])#18
Ext= np.array([0.501419, 0.635886, 0.770353, 0.904328, 1.03879, 1.22252, 1.55745,2.0239, 2.40661, 2.70067, 2.76125, 2.77061, 2.77997, 2.78933, 2.79869, 2.80805, 2.8174, 2.82676])##18

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
######################################################### 
pixel=float(45.0/512.0)
   
## 2022  1292
sou= [307.759 ,  191.924]  
lens=[305.916  , 185.771]   
msou= 15.867;  
mlens=17.868 
print "2022:  Lens_source_distance:  ",  np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
ds1=np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
print  "Delta_m1:  ",     abs(msou-mlens)#,  "Delta_m2:  ",  abs(mlens-mcomp)
input ("Enter a number ")    
print  "*************************************"

##2021
lens= [261.941 ,  245.252]##266.604,247.675] 
sou=[264.213 ,  251.164]  
msou= 9.888;  
mlens=11.817## 17.680    
print "2021:  Lens_source_distance:  ",  np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )
print  "Delta_m1:  ",     abs(msou-mlens)#,  "Delta_m2:  ",  abs(mlens-mcomp)
print "<v>:  ",  (np.sqrt((lens[0]-sou[0])**2.0  +( lens[1]-sou[1])**2.0 )- ds1) *45.0*1000.0/512.0

print "delta_dec(arcs):  ", abs(lens[0]-sou[0])*pixel,  "delta_ra(arcs):  ",  abs(lens[1]-sou[1])*pixel

input ("Enter a number ")    
    
    
    
###################Constant  ######################################
year=float(365.2422)
g=0.086
dx= 0.289

parcs= 30856775814913673.0
au= float( 1.495978707*(10.0**11.0))
G=6.67384*pow(10.,-11.0)#;// in [m^3/s^2*kg].
velocity= 3.0*pow(10.0,8.0)#;//velosity of light
Msun=1.98892*pow(10.,30)#; //in [kg].
phii=float(23.5*np.pi/180.0)
aup= float(1.0/206264.806247096);# AU/parsec
rad_arcs= float(180.0*3600.0/np.pi)## radian to arcs
omegae=float(2.0*np.pi/year);## earth angular velocity
cons= np.sqrt(4.0*G*Msun/(velocity * velocity * parcs * 1000.0) )
ti=np.zeros((2)); ex=np.zeros((2));  ey=np.zeros((2));  ez=np.zeros((2))  
########### OGLE-2011-BLG-1292 #####################################
fwhm= 6.6
Ftot=1744151.75## 5356901.50 # 7761073.00
dix= float(261.941-264.213)##pixel  North-South
diy= float(-1.0)*(245.252-251.164)##pixel  East-West 
ti[0]=float(159.0*24.0*3600.0+ 2.0*3600.0 + 6.0*60.0 + 24.0) ; ## t0 time with respect to Equinox 2011
ti[1]=float(3346.0*24.0*3600.0+ 3.0*3600.0 + 7.0*60.0 + 11.0);## 2021 image time  with .. ... 
parallax= 1.8813699680226963
sig_par= 0.6454810999999999
ms= 17.062 # I-band apparement magnitude of source 
sigms= 0.014     
dif= 1.929 
sigdif= 0.080
tE0= 5.614; 
sig_tE=0.039;
mus_x, mus_y= -9.437, -6.938 ## Gaia data
sig_musx, sig_musy= 0.470,  0.594

##########################################################################

mus= np.sqrt(mus_x*mus_x + mus_y*mus_y)
sig_mus= np.sqrt(sig_musx**2.0 + sig_musy**2.0)

Ds0=1.0/parallax 
sig_Ds= abs(sig_par/parallax/parallax)## kpc

ml= ms+dif;
sigml= np.sqrt( sigdif*sigdif + sigms*sigms)

distp=np.sqrt(dix*dix + diy*diy )## 13.156 pixel
dist=distp*pixel*1000.0 #mili arcs

sigma= np.sqrt(fwhm*fwhm/(g*Ftot * 2.0*np.log(2.0)) + dx*dx *np.log(2.0)/(2.0*np.pi* fwhm*fwhm))##pixel
sig_dp= sigma*np.sqrt(2.0)#pixel
sig_d= sig_dp*pixel*1000.0 #mili arcs

tim = float((ti[1]-ti[0])/3600.0/24.0/year)## total time in year
sig_t=float(1.0/year) ## error in time in day

mu= float(dist/tim)### mili arcs/year
sig_mu= mu* np.sqrt((sig_dp/distp)**2.0 + (sig_t/tim)**2.0) 
mux=float(dix*pixel*1000.0/tim)##marcs/year
muy=float(diy*pixel*1000.0/tim)##marcs/year
sig_mux= abs(mux)* np.sqrt( (sigma/dix)**2.0  +  (sig_t/tim)**2.0) 
sig_muy= abs(muy)* np.sqrt( (sigma/diy)**2.0  +  (sig_t/tim)**2.0) 


##########################################################################


ns= 300
Ml=np.zeros((ns));  Dl= np.zeros((ns))
for i in range(ns): 
    Ml[i]=pow(10.0,-1.11 + i*1.11/ns ) ##float(0.05 + (1.0-0.05)*i/ns) 
    Dl[i]=float(0.0005 +(Ds0-0.001)*i/ns) ### [0.1,7.9]  
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
###############################################################
if(len(dl1)<len(dl4)): num=len(dl1);
else :                 num=len(dl4)        
tEE=np.zeros((num,4))            
tEE[:, 0]=dl1[:num] 
tEE[:, 1]=ml1[:num]    
tEE[:, 2]=dl4[:num]           
tEE[:, 3]=ml4[:num]    
fil2=open('./file1292.dat',"w"); 
np.savetxt(fil2,tEE.reshape((-1,4)) ,fmt="%.6f     %.6f    %.6f      %.6f")  
fil2.close()  



############################################################################       
dd=0.0;  mm=0.0; count=0 
for i in range(len(dl1)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl1[i]-dl4[j])**2.0 + (ml1[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl1[i];   mm+=ml1[i];  count+=1.0 
DL= float(dd/count)  
ML= float(mm/count)                    
print "<DL>:   ",  DL,  "  <ML>:    ",  ML          
############################################# 
dd=0.0;  mm=0.0; count=0 
for i in range(len(dl2)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl2[i]-dl4[j])**2.0 + (ml2[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl2[i];   mm+=ml2[i];  count+=1.0 
DL1= float(dd/count)  
ML1= float(mm/count)                    
print "<DL>_Down:   ",  DL1,  "  <ML>_Down:    ",  ML1         
############################################# 
dd=0.0;  mm=0.0; count=0 
for i in range(len(dl3)):  
    for j in range(len(dl4)):  
       dis=np.sqrt((dl3[i]-dl4[j])**2.0 + (ml3[i]-ml4[j])**2.0)
       if(dis<0.01):
           dd+= dl3[i];   mm+=ml3[i];  count+=1.0 
DL2= float(dd/count)  
ML2= float(mm/count)                    
print "<DL>_Upper:   ",  DL2,  "  <ML>_Upper:    ",  ML2      
           
           



 
############################################################################       

       
print "******************  OGLE-2012-BLG-1292 **********************"
print "Ds:  ", Ds0 ,  " sig_Ds:  ", sig_Ds
print "sigDL:  ", DL1,  DL2 ,    " sihmass:  ", DL1, DL2
print "mu(marcs/yrs):  ", mu, "sig_mu: ",  sig_mu 
print  "M_app(lens)[mag]:    ",  ml , "Error_Map_lens:  ",  sigml 
print "tE(days),  ", tE0, "   sig_tE: ", sig_tE
print "Dl(kpc):  ",  DL,  "  Mass(Solar mass):  ",  ML           
print "tetE (miliarcs):  ",  tetE
print "**************************************************************"


###########################################################################
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
plt.xlim([0.05,Ds0])
plt.ylim([0.08,np.max(Ml)])
plt.title(r"OGLE-2012-BLG-1292", fontsize=18)
plt.text(0.07,0.8,r"$t_{\rm E}(\rm{days})=$"+str(round(tE0,2))+ r"$\pm$" +str(round(sig_tE,2)), fontsize=19, color="r")
plt.text(0.07,0.7,r"$m_{\rm I,~\rm{lens}}(\rm{mag})=$"+str(round(ml,2))+r"$\pm$"+str(round(sigml,2)),fontsize=19, color="b")
fig=plt.gcf()
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig.savefig("1292eventc.jpg",dpi=200)
py.clf()
##############################################################################    
        
xrel=float(DL/Ds0)
tetE=cons*np.sqrt(ML*(1.0-xrel)/xrel/Ds0) * 180.0/np.pi * 3600.0 *1000.0 ## mili arcs
sigtetE=tetE*(sig_tE/tE0 + sig_mu/mu)  
sig_Dl=np.sqrt((DL1-DL)**2.0 +  (DL2-DL)**2.0 ) ## 1.86896## kpc
pi_rel= float(1.0/DL - 1.0/Ds0)## mili arcs
sig_pi= np.sqrt((sig_Dl/DL/DL)**2.0 +  (sig_Ds/Ds0/Ds0)**2.0 )## milli arcs

  


       
       
  
for i in range(2):  
    ex[i]= aup*pi_rel*rad_arcs* np.cos(omegae*ti[i]) ; ##projection in the line of sight  [marcs]
    ey[i]= aup*pi_rel*rad_arcs* np.sin(omegae*ti[i])*np.cos(phii)#projected in RA. direction [marcs]
    ez[i]= aup*pi_rel*rad_arcs* np.sin(omegae*ti[i])*np.sin(phii)#projected in DEC direction [marcs]
del_earx= float(ez[1]-ez[0]) ### marcs
del_eary= float(ey[1]-ey[0]) ### marcs   

Del_helx= float(dix*pixel*1000.0 - del_earx) ## marcs
Del_hely= float(diy*pixel*1000.0 - del_eary) ## marcs
sig_hel=np.sqrt( (sigma*pixel*1000.0)**2.0 +  sig_pi**2.0)## marcs

muxh=float(Del_helx/tim)##marcs/year
muyh=float(Del_hely/tim)##marcs/year
sig_muxh= abs(muxh)* np.sqrt( (sig_hel/Del_helx)**2.0  +  (sig_t/tim)**2.0) 
sig_muyh= abs(muyh)* np.sqrt( (sig_hel/Del_hely)**2.0  +  (sig_t/tim)**2.0)      
     
mulx= mus_x + muxh 
muly= mus_y + muyh

sig_mulx= np.sqrt( sig_muxh**2.0  +  sig_musx**2.0 ) 
sig_muly= np.sqrt( sig_muyh**2.0  +  sig_musy**2.0 )      
   
  
print "******************************** OGLE-2011-BLG-1292 ****************"
print "Distance_x (marcs): ",   float(dix*pixel*1000.0), "\pm" , sigma*pixel*1000.0
print "Distance_Y (marcs):  ",  float(diy*pixel*1000.0), "\pm" , sigma*pixel*1000.0
print "dis_size +-  sigma (pixel):  ",    distp ,    sig_dp 
print "dis_size +-  sigma (marcs):  ",    dist ,     sig_d 
print "time     +-  sigma (year):  ",     tim ,      sig_t
print  "mu      +-  dmu: (marcs/year):",  mu ,      sig_mu
print  "mu_N(geo)     +-  dmu_N: (marcs/year):",  mux ,      sig_mux
print  "mu_E(geo)     +-  dmu_E: (marcs/year):",  muy ,      sig_muy
print "\nDisx_helo:(marcs):  ",  Del_helx,  "\pm:  ",  sig_hel
print "Disx_helo:(marcs):  ",  Del_hely,  "\pm:  ",  sig_hel
print  "mu_N(helo)     +-  dmu_N: (marcs/year):",   muxh,      sig_muxh
print  "mu_E(helo)     +-  dmu_E: (marcs/year):",   muyh,      sig_muyh

print  "\nmu_NLens(helo)     +-  dmu_N: (marcs/year):",  mulx,    sig_mulx
print  "mu_ELens(helo)     +-  dmu_E: (marcs/year):",    muly,    sig_muly


mull=np.sqrt(mulx*mulx + muly*muly)
errm= np.sqrt(sig_mulx*sig_mulx +  sig_muly*sig_muly)
print "mu_lens(size):  ",  mull ,    errm
print "V_lens(km/s): ",    mull * DL *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_up: ",    (mull) * DL1 *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_down: ",    (mull) * DL2 *au*0.001/(year*24.0*3600.0)

print "\n\nmag_lens(I-band) \pm sigma:  ",   ml ,    sigml 
 
print "\n\nEarth displacement during 9.94 years:  RA_projected[marcs]   ",  del_earx
print "   Earth displacement during 9.94 years:  DEC_projected[marcs]:   " , del_eary
print "Relative parallax (mili arcs):  ", float(pi_rel* aup* rad_arcs),  " Error:  ",  sig_pi
print "Extinction at the lens position:    ",  extinction(DL)
print "**********************************************************************"


param= np.array([tE0,sig_tE , mu, sig_mu , ml,sigml, DL, DL1, DL2 , ML, ML1, ML2,  tetE, sigtetE, mus, sig_mus ,Ds0, sig_Ds, mull, errm])
filf=open('./Tab_paper1292.dat',"w"); 
np.savetxt(filf,param.reshape((-1,20)), fmt="$OGLE1292$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f \pm %.2f$ & $%.2f\pm%.2f \pm %.2f$ & $%.2f\pm%.2f$ & $%.2f\pm%.2f$  & $%.2f\pm%.2f$ & $%.2f\pm%.2f$")
filf.close()       




import numpy as np
import matplotlib.pyplot as plt
import pylab as py
from matplotlib import rcParams
from numpy import ma
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 11.5
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
import math
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib as mpl
from matplotlib import gridspec
from matplotlib.ticker import FormatStrFormatter

###=============================================================================
'''
nn=1526
data=np.zeros(( nn , 7 ))
data=np.loadtxt('./OGLE.txt')
tes=0
for i in range(nn):  
    if(data[i,0]<10.0):   
        tes+=1.0
print "No.  OGLE events with tE<10:  ",  tes
print "<tE:  >",  np.mean(data[:,0]),  np.std(data[:,0])
'''

# Equinox https://data.giss.nasa.gov/cgi-bin/ar5/srevents.cgi 
# time calcu:  https://www.timeanddate.com/date/durationresult.html?d1=20&m1=3&y1=2011&d2=17&m2=5&y2=2021&h1=23&i1=31&s1=0&h2=6&i2=49&s2=47

year=float(365.2422)
g=0.086
dx= 0.289
pixel=float(45.0/512.0)
omegae=float(2.0*np.pi/year);## earth angular velocity
au=   float( 1.495978707*(10.0**11.0))
phii=float(23.5*np.pi/180.0)
aup= float(1.0/206264.806247096);# AU/parsec
rad_arcs= float(180.0*3600.0/np.pi)## radian to arcs
parcs= 30856775814913673.0
ti=np.zeros((2)); ex=np.zeros((2));  ey=np.zeros((2));  ez=np.zeros((2))  


####################################  OGLE-2011-BLG-0433 #####################################
fwhm= 5.5
Ftot=279851.875## 7761073.00  ????
dix= float(276.504-266.604) ##pixel  North-Sourth (declination)
diy= float(-1.0)*(239.003-247.675)##pixel East-West (right ascention)
ti[0]=float(75.0*24.0*3600.0+ 18.0*3600.0 + 6.0*60.0 + 12.0) ; ## t0 time with respect to Exoinox 2011
ti[1]=float(3710.0*24.0*3600.0+ 7.0*3600.0 + 18.0*60.0 +47.0);## 2021 image time  with .. ... 

mus_x, mus_y= -8.189, -1.085 ## Gaia data
sig_musx, sig_musy= 0.302,  0.422

Dl= 0.38093409367183795  #0.134288492477##2.2446932352#k parcs
Dl1=0.37606427326097047 #0.125274025242 ##  1.04818889308 
Dl2=0.38729562341312296 #0.142122037745### 2.96741412331
Ds=8.0##  0.424498929378  #k parcs
sig_Dl=np.sqrt((Dl1-Dl)**2.0 +  (Dl2-Dl)**2.0 ) ## 1.86896## kpc
sig_Ds=0.5## 0.0681008448672


ms=    15.798 # I-band apparement magnitude of source 
sigms= 0.000    
dif=   2.984 
sigdif=0.037 


distp=np.sqrt(dix*dix + diy*diy )#pixel
sigma= np.sqrt( fwhm*fwhm/(g*Ftot * 2.0*np.log(2.0)) + dx*dx *np.log(2.0)/(2.0*np.pi* fwhm*fwhm))##pixel
sig_dp= sigma*np.sqrt(2.0)#pixel
dist=distp*pixel*1000.0#mili arcs
sig_d= sig_dp*pixel*1000.0#mili arcs
tim = float((ti[1]-ti[0])/3600.0/24.0/year)## total time in year
sig_t= float(1.0/year) ## error in time in day
mu= float(dist/tim)### mili arcs/year
sig_mu= abs(mu)* np.sqrt((sig_dp/distp)**2.0 + (sig_t/tim)**2.0) 
mux=float(dix*pixel*1000.0/tim)##marcs/year
muy=float(diy*pixel*1000.0/tim)##marcs/year
sig_mux= abs(mux)* np.sqrt( (sigma/dix)**2.0  +  (sig_t/tim)**2.0) 
sig_muy= abs(muy)* np.sqrt( (sigma/diy)**2.0  +  (sig_t/tim)**2.0) 
pi_rel= float(1.0/Dl - 1.0/Ds)## mili arcs
sig_pi= np.sqrt((sig_Dl/Dl/Dl)**2.0 +  (sig_Ds/Ds/Ds)**2.0 )## milli arcs
ml= ms+dif;
sigml= np.sqrt( sigdif*sigdif + sigms*sigms )  
    
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
   
   
print ("******************************** OGLE-2011-BLG-0433 ****************")
print ("Distance_x (marcs): ",   float(dix*pixel*1000.0), "\pm" , sigma*pixel*1000.0)
print ("Distance_Y (marcs):  ",  float(diy*pixel*1000.0), "\pm" , sigma*pixel*1000.0)
print ("dis_size +-  sigma (pixel):  ",    distp ,    sig_dp )
print ("dis_size +-  sigma (marcs):  ",    dist ,     sig_d )
print ("time     +-  sigma (year):  ",     tim ,      sig_t)
print  ("mu      +-  dmu: (marcs/year):",  mu ,      sig_mu)
print  ("mu_N(geo)     +-  dmu_N: (marcs/year):",  mux ,      sig_mux)
print  ("mu_E(geo)     +-  dmu_E: (marcs/year):",  muy ,      sig_muy)
print ("\nDisx_helo:(marcs):  ",  Del_helx,  "\pm:  ",  sig_hel)
print ("Disx_helo:(marcs):  ",  Del_hely,  "\pm:  ",  sig_hel)
print  ("mu_N(helo)     +-  dmu_N: (marcs/year):",   muxh,      sig_muxh)
print  ("mu_E(helo)     +-  dmu_E: (marcs/year):",   muyh,      sig_muyh)
print  ("\nmu_NLens(helo)     +-  dmu_N: (marcs/year):",  mulx,    sig_mulx)
print  ("mu_ELens(helo)     +-  dmu_E: (marcs/year):",    muly,    sig_muly)



mull=np.sqrt(mulx*mulx + muly*muly)
errm= np.sqrt(sig_mulx*sig_mulx +  sig_muly*sig_muly)
print ("mu_lens(size):  ",  mull ,    errm)
print ("V_lens(km/s): ",    mull * Dl *au*0.001/(year*24.0*3600.0) )
print ("V_lens(km/s)_up: ",    (mull) * Dl1 *au*0.001/(year*24.0*3600.0) )
print ("V_lens(km/s)_down: ",    (mull) * Dl2 *au*0.001/(year*24.0*3600.0))



print ("\n\nmag_lens(I-band) \pm sigma:  ",   ml ,    sigml )
print ("\n\nEarth displacement during 9.94 years:  RA_projected[marcs]   ",  del_earx)
print ("   Earth displacement during 9.94 years:  DEC_projected[marcs]:   " , del_eary)
print ("\n\nEarth displacement during 9.94 years:  RA_projected[au]   ",  del_earx/pi_rel)
print ("   Earth displacement during 9.94 years:  DEC_projected[au]:   " , del_eary/pi_rel)
print ("Relative parallax (mili arcs):  ", float(pi_rel* aup* rad_arcs),  " Error:  ",  sig_pi)
print ("**********************************************************************")

input("Enter a number !!!")



###########################  OGLE-2012-BLG-1276 ######################################
fwhm= 5.0
Ftot=5356901.50 # 7761073.00
dix= float(243.631 - 248.979)##pixel North-South
diy= float(-1.0)*(242.506 - 244.129)##pixel East-West
ti[0]=float(156.0*24.0*3600.0+ 22.0*3600.0 + 30.0*60.0 + 24.0) ; ## t0 time with respect to Exoinox 2011
ti[1]=float(3352.0*24.0*3600.0+ 1.0*3600.0 + 51.0*60.0 + 24.0);## 2021 image time  with .. ... 


mus_x, mus_y= -0.895, -4.892 ## Gaia data
sig_musx, sig_musy= 0.106,  0.158

Dl=2.29429074154##2.2446932352#k parcs
Dl1=1.16924343542##  1.04818889308 
Dl2=2.99115394043### 2.96741412331
Ds=  4.87725948097 #k parcs
sig_Dl=np.sqrt((Dl1-Dl)**2.0 +  (Dl2-Dl)**2.0 ) ## 1.86896## kpc
sig_Ds=  2.78666252633

ms= 16.148 # I-band apparement magnitude of source 
sigms= 0.025     
dif= 3.736 
sigdif= 0.079


distp=np.sqrt(dix*dix + diy*diy )#pixel
sigma= np.sqrt( fwhm*fwhm/(g*Ftot * 2.0*np.log(2.0)) + dx*dx *np.log(2.0)/(2.0*np.pi* fwhm*fwhm))##pixel
sig_dp= sigma*np.sqrt(2.0)#pixel
dist=distp*pixel*1000.0#mili arcs
sig_d= sig_dp*pixel*1000.0#mili arcs
tim = float((ti[1]-ti[0])/3600.0/24.0/year)## total time in year
sig_t= float(1.0/year) ## error in time in day
mu= float(dist/tim)### mili arcs/year
sig_mu= abs(mu)* np.sqrt((sig_dp/distp)**2.0 + (sig_t/tim)**2.0) 
mux=float(dix*pixel*1000.0/tim)##marcs/year
muy=float(diy*pixel*1000.0/tim)##marcs/year
sig_mux= abs(mux)* np.sqrt( (sigma/dix)**2.0  +  (sig_t/tim)**2.0) 
sig_muy= abs(muy)* np.sqrt( (sigma/diy)**2.0  +  (sig_t/tim)**2.0) 
pi_rel= float(1.0/Dl - 1.0/Ds)## mili arcs
sig_pi= np.sqrt((sig_Dl/Dl/Dl)**2.0 +  (sig_Ds/Ds/Ds)**2.0 )## milli arcs
ml= ms+dif;
sigml= np.sqrt( sigdif*sigdif + sigms*sigms )  
    
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
   
  
print "******************************** OGLE-2011-BLG-1276 ****************"
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
print "V_lens(km/s): ",    mull * Dl *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_up: ",    (mull) * Dl1 *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_down: ",    (mull) * Dl2 *au*0.001/(year*24.0*3600.0)



print "\n\nmag_lens(I-band) \pm sigma:  ",   ml ,    sigml 
 
print "\n\nEarth displacement during 9.94 years:  RA_projected[marcs]   ",  del_earx
print "   Earth displacement during 9.94 years:  DEC_projected[marcs]:   " , del_eary
print "Relative parallax (mili arcs):  ", float(pi_rel* aup* rad_arcs),  " Error:  ",  sig_pi
print "**********************************************************************"

input("Enter a number !!!")


##################################  OGLE-2012-BLG-1292 #########################################
fwhm= 6.6
Ftot=1744151.75## 5356901.50 # 7761073.00
dix= float(261.941-264.213)##pixel  North-South
diy= float(-1.0)*(245.252-251.164)##pixel  East-West 
ti[0]=float(159.0*24.0*3600.0+ 2.0*3600.0 + 6.0*60.0 + 24.0) ; ## t0 time with respect to Exoinox 2011
ti[1]=float(3346.0*24.0*3600.0+ 3.0*3600.0 + 7.0*60.0 + 11.0);## 2021 image time  with .. ... 

mus_x, mus_y= -9.437, -6.938 ## Gaia data
sig_musx, sig_musy= 0.470,  0.594


### CHANGED
Dl= 0.389735995869##  0.378881403874#kpc
Ds=  0.53152756608 #kpc   
Dl1=0.264977400237## 0.258826096481
Dl2=0.515895576035### 0.498013858069

sig_Dl=np.sqrt((Dl1-Dl)**2.0 +  (Dl2-Dl)**2.0 ) ## 1.86896## kpc
sig_Ds=0.182362323129#kpc

ms= 17.062 # I-band apparement magnitude of source 
sigms= 0.014     
dif= 1.929 
sigdif=0.080




distp=np.sqrt(dix*dix + diy*diy )#pixel
sigma= np.sqrt( fwhm*fwhm/(g*Ftot * 2.0*np.log(2.0)) + dx*dx *np.log(2.0)/(2.0*np.pi* fwhm*fwhm))##pixel
sig_dp= sigma*np.sqrt(2.0)#pixel
dist=distp*pixel*1000.0#mili arcs
sig_d= sig_dp*pixel*1000.0#mili arcs
tim = float((ti[1]-ti[0])/3600.0/24.0/year)## total time in year
sig_t= float(1.0/year) ## error in time in day
mu= float(dist/tim)### mili arcs/year
sig_mu= abs(mu)* np.sqrt((sig_dp/distp)**2.0 + (sig_t/tim)**2.0) 
mux=float(dix*pixel*1000.0/tim)##marcs/year
muy=float(diy*pixel*1000.0/tim)##marcs/year
sig_mux= abs(mux)* np.sqrt( (sigma/dix)**2.0  +  (sig_t/tim)**2.0) 
sig_muy= abs(muy)* np.sqrt( (sigma/diy)**2.0  +  (sig_t/tim)**2.0) 
pi_rel= float(1.0/Dl - 1.0/Ds)## mili arcs
sig_pi= np.sqrt((sig_Dl/Dl/Dl)**2.0 +  (sig_Ds/Ds/Ds)**2.0 )## milli arcs
ml= ms+dif;
sigml= np.sqrt( sigdif*sigdif + sigms*sigms )  
    
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
print "V_lens(km/s): ",    mull * Dl *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_up: ",    (mull) * Dl1 *au*0.001/(year*24.0*3600.0)
print "V_lens(km/s)_down: ",    (mull) * Dl2 *au*0.001/(year*24.0*3600.0)




print "\n\nmag_lens(I-band) \pm sigma:  ",   ml ,    sigml 
 
print "\n\nEarth displacement during 9.94 years:  RA_projected[marcs]   ",  del_earx
print "   Earth displacement during 9.94 years:  DEC_projected[marcs]:   " , del_eary
print "Relative parallax (mili arcs):  ", float(pi_rel* aup* rad_arcs),  " Error:  ",  sig_pi
print "**********************************************************************"

input("Enter a number !!!")












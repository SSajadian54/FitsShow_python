#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include "VBBinaryLensingLibrary.h"
using namespace std;

const int YZ=3615;
const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const double KP=3.08568025*pow(10.,19); // in meter.
const double G= 6.67384*pow(10.,-11.0);// in [m^3/s^2*kg].
const double velocity= 3.0*pow(10.0,8.0);//velosity of light
const double M_sun=1.98892*pow(10.,30); //in [kg].
const double Rsun= 6.957*pow(10.0,8.0); ///solar radius [meter]
const double Mjupiter=1.898*pow(10,27); 
const double logT_sun=log10(5778.0);
const double vro_sun=226.0;
const double AU=1.4960*pow(10.0,11.0);
const double year=365.24222; 
const double binary_fraction=double(2.0/3.0);

const int Nef= 151;
const int Nfi= 144; 
const int Nog= 980;
///============================ Besancon constant ==============///
const double R_sun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]=  {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]=  {3.1,2.5,3.1,3.1};
const int N1=19744,N2=43470, N3=4404, N4=2597;///CMD_BESANCON, thinD, bulge, thickD, halo
///============================== MOA observation constants ====================
const int M=4; ///filters UBVI
const double wav[M]=  {0.4353,0.5477,0.6349,0.8797};///in[micrometer] BVRI
const double sigma[M]={0.022,0.022,0.02,0.025};//MOAاستفاده از مقاله کاردلی 
const double FWHM[M]= { 3.5*0.26, 3.5*0.26, 3.5*0.26, 3.5*0.26};//MOA[arcsec] 5.5,5.0,4.5,4.0 times Pixel size (0.58'')////BVRI
const double AlAv[M]= {1.555,1.332,1.009,0.600};///From besancon model[UBVI]
const double Avks=double(8.20922); 
const double VSunR =11.1;
const double VSunT =vro_sun*(1.00762+0.00712)+ 12.24;
const double VSunZ =7.25;
///======================================================
struct source{
    int nums,struc, cl;
    double Ds,TET,FI;
    double lon,lat,col, type, mass;
    double logl, logg,Teff, Rs, ro_star, magni;
    double od_disk,od_ThD,od_bulge,od_halo,opd;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs,nstart;
    double Mab[M],Map[M];
    double Fluxb[M],magb[M]; 
    double blend[M],nsbl[M]; 
    double SV_n1, LV_n1, VSun_n1;
    double SV_n2, LV_n2, VSun_n2;
    double mus1, mus2, mul1, mul2;
    double ext[M], age, deltao;
};
struct lens{
    int numl,struc, cl;
    double mul, Rl; 
    double Ml,Dl, vl, vs, Vt, xls;
    double rhomaxl, tE, RE;
    double type, Teff, logl, logg, age, u0;
    double Mapl[M]; 
    double pt1, pt2, dt;
};
struct CMD{
    double Teff_d[N1],logl_d[N1],Mab_d[M][N1],mass_d[N1], type_d[N1], metal_d[N1],age_d[N1], gra_d[N1],Rs_d[N1]; int cl_d[N1];//1disk
    double Teff_b[N2],logl_b[N2],Mab_b[M][N2],mass_b[N2], type_b[N2], metal_b[N2],age_b[N2], gra_b[N2],Rs_b[N2]; int cl_b[N2];//bulge
    double Teff_t[N3],logl_t[N3],Mab_t[M][N3],mass_t[N3], type_t[N3], metal_t[N3],age_t[N3], gra_t[N3],Rs_t[N3]; int cl_t[N3];//2disk
    double Teff_h[N4],logl_h[N4],Mab_h[M][N4],mass_h[N4], type_h[N4], metal_h[N4],age_h[N4], gra_h[N4],Rs_h[N4]; int cl_h[N4];//halo
};
struct extinc{
   double dis[100];///distance 
   double Extks[100];///ks-band extinction
   double Ai[M],Av,Aks,AI;
   double exks;
};
///===================== FUNCTION ==============================================
void read_cmd(CMD & cm);
void func_source(source & s,CMD & cm,extinc & ex);
void func_lens(lens & l, source & s, CMD & cm, extinc & ex);
void vrel(source & s,lens & l);
void Disk_model(source & s);
void optical_depth(source & s);
int Extinction(extinc & ex,source & s);
double Interpol(double ds, extinc & ex);
double ErrorCal(double);
double RandN(double , double);
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    time_t _timeNow;
    unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///==============================================================//
///                                                              //
///                  Main program                                //
///                                                              //
///==============================================================//
int main()
{
//================================================
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
///=================================================

     CMD cm;
     source s;
     lens l;
     extinc ex;
     read_cmd(cm);
     
     
     
     VBBinaryLensing vbb;
     printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
     vbb.Tol=1.e-4;
     vbb.LoadESPLTable("./files/ESPL.tbl");
     vbb.a1=0.0;  
    
    
    
    FILE* fil1;  FILE* fil2;  FILE* fil3;  FILE* ogle; //FILE *efimoa;
    char filename1[40], filename2[40];
   // efimoa=fopen("./files/effi_black_MOA.txt","r"); 
   // double tef[Nef]={0.0}; double efi[Nef]={0.0};
    //for(int i=0; i<Nef; ++i){
    //fscanf(efimoa,"%lf    %lf \n",&tef[i],&efi[i]);}
    //fclose(efimoa);
    fil1=fopen("./files/distribution/dist_SDML.txt","w");
    ogle=fopen("./files/OGLE_delta.dat","r");
    double dto[Nog], dmo, magg;  
    for(int i=0; i<Nog; ++i){fscanf(ogle,"%lf   %lf   %lf \n", &dto[i], &dmo, &magg);}
    fclose(ogle);
   
   
   
    int  ndt=0, flagm=0, nsimu=0, ndet;  
    double Astar, sigA, tt, maga, Efi2, u; 
    double tim, timp; 
    double chi0=0.0, chi1=0.0, dchi, thre, night;
    double flag[3000]; 
    
    double cade=double(20.0/60.0/24.0);
    
   
    s.lat=-3.5;
    double lonn=7.0;
    cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
    if(lonn<=0.0)   s.lon=360.0+lonn;
    else            s.lon=lonn;
    s.TET=(360.0-s.lon)/RA;///radian
    s.FI=s.lat/RA;///radian
    if(Extinction(ex,s)==1){
    Disk_model(s);
    ndet=0;  nsimu=0; 


    for(nsimu=0; nsimu<10000000; ++nsimu){
    

    do{
    func_source(s, cm, ex);
    optical_depth(s);
    func_lens(l, s, cm, ex);
    }while(l.tE<0.001  or l.tE>12.0 or s.type>=8.0 or s.cl>5 or l.Vt<10.0 or l.Vt>350.0); 
    
    
    
    Efi2=double(((double)rand()/(double)(RAND_MAX+1.))*1.0);
    Astar=vbb.ESPLMag2(l.u0 , s.ro_star);
    maga=double(s.magb[3] -2.5*log10(s.blend[3]*Astar + 1.0- s.blend[3]));
    if(maga>12.0 and s.magb[3]<21.0 and Efi2<double(s.blend[3])){
    
    
    flagm=0; 
    if(nsimu%10000==0){
    flagm=1;
    sprintf(filename1,"./files/distribution/%c%c%d.dat",'d','_',int(nsimu) );
    fil2=fopen(filename1,"w");
    sprintf(filename2,"./files/distribution/%c%c%d.dat",'l','_',int(nsimu) );
    fil3=fopen(filename2,"w");}
   
  
    chi0=0.0, chi1=0.0, thre=0.0; ndt=0; 
    timp=0;  
    night= double((double)rand()*(Nog-20)/(double)(RAND_MAX+1.))*24.0;
    
    for(double tim=l.pt1;  tim<l.pt2;  tim=tim+l.dt){
    timp +=l.dt; 
    night+=double(l.dt*24.0); 
    if(night>24.0)  night=double(night-24.0); 
    
    u= sqrt((tim/l.tE)*(tim/l.tE) + l.u0*l.u0 ); 
    Astar=vbb.ESPLMag2(u , s.ro_star);
    Astar=fabs(s.blend[3]*Astar+1.0-s.blend[3]); 
    maga= s.magb[3]-2.5*log10(Astar); 
    if(flagm==1) fprintf(fil3, "%.5lf    %.5lf\n", tim, Astar); 
    
    if(timp>cade){
    timp=timp-cade;

    if(maga>12.0 and maga<21.0 and night<6.5){
    dmo=ErrorCal(maga); 
    sigA=fabs(pow(10.0,-0.4*dmo)-1.0)*Astar;
    tt= RandN(1.0, 2.5); 
    chi0 +=fabs( (sigA*tt)*(sigA*tt)/(sigA*sigA) );  
    chi1 +=fabs( (Astar+sigA*tt - 1.0)*(Astar+sigA*tt - 1.0 )/(sigA*sigA) );
    
    flag[ndt]=-1.0; 
    if(double(Astar+sigA*tt-1.0)>fabs(4.0*sigA) )    flag[ndt]=1.0; 
    if(ndt>2 and double(flag[ndt]+flag[ndt-1]+flag[ndt-2]+flag[ndt-3])>3.0)  thre=1.0; 
    if(flagm==1) fprintf(fil2,"%.5lf  %.5lf  %.5lf\n", tim, Astar+sigA*tt, sigA );
    ndt+=1; }
  
    if(dmo<0.0 or sigA<0.0 or maga<0.0 or Astar<1.0 or ndt<0 or chi0<0.0 or chi1<0.0 or s.magb[3]<12.0 or ndt>2998){
    cout<<"ERROR !!! dmo:  "<<dmo<<"\t sigA:  "<<sigA<<"\t magg:  "<<magg<<endl;
    cout<<"magb[3]:    "<<s.magb[3]<<"\t Astar:  "<<Astar<<"\t blending:  "<<s.blend[3]<<"\t u0:  "<<l.u0<<endl;
    cout<<"Astar:   "<<Astar<<"\t ndt:  "<<ndt<<"\t magb:  "<<s.magb[3]<<endl; int yye; cin>>yye;}
    }}
   
    if(flagm==1){ fclose(fil2); fclose(fil3); }
    dchi=fabs(chi0-chi1);
    if(dchi>800.0 and thre>0.5){
    ndet+=1; 
    fprintf(fil1,"%d %.2lf %.9lf  %.3lf  %.3lf  %.5lf  %.3lf  %.3lf  %.4lf  %.4lf  %.7lf  %.6lf  %.1lf  %.1lf %d\n",
    nsimu,l.type,l.Ml,l.Dl,s.Ds,l.tE,l.Vt,l.u0,s.blend[3], s.magb[3],s.ro_star,l.RE/AU, thre, dchi, ndt);
   
    
   if(fabs(double(s.Map[0]-s.Mab[0]-s.ext[0])-double(s.Map[3]-s.Mab[3]-s.ext[3]))>0.1 or s.mass<0.0 or l.tE<0.0 or l.vl<0.0 or l.Vt<0.0 or l.Dl<0.0 or s.Ds<0.0 or s.Ds>20.0 or l.Dl>s.Ds or fabs(l.u0)>1.5 or l.Ml<0.0 or s.ro_star<0.0 or s.Rs<=0.0 or s.type>8.0 or s.type<2.0 or (s.cl==5 and int(s.type)>7) or (s.cl==6 and int(s.type)!=9) or (s.cl<5 and int(s.type)==9)){
   cout<<"ERROR distance medule(1): "<<s.Map[0]-s.Mab[0]-s.ext[0]<<"\t distance medule(3): "<<s.Map[3]-s.Mab[3]-s.ext[3]<<endl;
   cout<<"Map[0]:  "<<s.Map[0]<<"\t Mab[0]:  "<<s.Mab[0]<<"\t extinc:  "<<s.ext[0]<<endl;
   cout<<"Map[3]:  "<<s.Map[3]<<"\t Mab[3]:  "<<s.Mab[3]<<"\t extinc:  "<<s.ext[3]<<endl;
   cout<<"mass: "<<s.mass<<"\t tE: "<<l.tE<<"\t vl: "<<l.vl<<"\t Vt: "<<l.Vt<<endl; 
   cout<<"Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<"\t u0: "<<l.u0<<"\t Ml: "<<l.Ml<<endl; 
   cout<<"ro_star: "<<s.ro_star<<"\t Rs: "<<s.Rs<<endl;
   cout<<"type: "<<s.type<<"\t cl: "<<s.cl<<endl; int uue; cin>>uue;}}}
  
   cout<<"**********************************************"<<endl; 
   cout<<"nsimu:  "<<nsimu<<"\t ndet:  "<<ndet<<endl; 
   cout<<"tE:  "<<l.tE<<"\t Vt:  "<<l.Vt<<"\t Ml:  "<<l.Ml<<endl;
   cout<<"typ_lens:  "<<l.type<<"\t RE:  "<<l.RE/AU<<"\t ro_star:  "<<s.ro_star<<endl;
   cout<<"**********************************************"<<endl;  }
   }
    fclose(fil1); 
    time(&_timeNow);
    printf("END_TIME >>>>>>>>  %s ",ctime(&_timeNow));
    return(0);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double ErrorCal(double mag){
   double err; 

   if(mag<15.0 ) err=0.003;
   if(mag>=15.0) err=0.00862339038 + 0.00682867290*(mag-16.6494077) +  0.00227422710*(mag-16.6494077)*(mag-16.6494077);  
   if(mag<19 and err >0.095) {cout<<"Error(OGLE) error:  "<<err<<"\t mag:  "<<mag<<endl;  int uue; cin>>uue;  }
   
   return(err);
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex)
{
  double F=-1.0;
  if(ds<ex.dis[0])       F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm){
//// mass, teff, Age, logL, logg,  metal, R*, U, B, V, I,  K,  Cl, Type
///// 0     1     2     3     4    5      6   7  8  9  10  11  12   13 
    int yye, uui,h, k1, k2, g; double metal, mk;   
    char filename[40];
    FILE *fp2;  
////=================================== THIN DISK ==============================
    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c%c.dat",'C','M','D','T','i');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&cm.age_d[j],&cm.logl_d[j],&cm.gra_d[j],&cm.metal_d[j],&cm.Rs_d[j],
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[3][j],&mk,&cm.cl_d[j],&cm.type_d[j]);
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.gra_d[j]>6.0 or cm.Teff_d[j]<0.0 or cm.metal_d[j]>0.12 or cm.age_d[j]>10.0 or cm.age_d[j]<0.0 or
    cm.cl_d[j]>6 or cm.type_d[j]>=9.0 or cm.type_d[j]<2.0 or (cm.cl_d[j]==5 and int(cm.type_d[j])>7) or 
    (cm.cl_d[j]==6 and int(cm.type_d[j])!=9) or (cm.cl_d[j]<5 and int(cm.type_d[j])==9)){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; 
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;




////=================================== Galactic Bulge ===========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c2.dat",'C','M','D','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&cm.age_b[j],&cm.logl_b[j],&cm.gra_b[j],&cm.metal_b[j],&cm.Rs_b[j],
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[3][j],&mk,&cm.cl_b[j],&cm.type_b[j]);
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or cm.age_b[j]>10.0 or cm.metal_b[j]>0.15 or cm.cl_b[j]>5 or 
    cm.type_b[j]>8.0 or (cm.cl_b[j]==5 and int(cm.type_b[j])>7) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) or 
    (cm.cl_b[j]<5 and int(cm.type_b[j])==9)){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;



////==================================== THICK DISK ==========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c%c.dat",'C','M','D','T','k');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&cm.age_t[j],&cm.logl_t[j],&cm.gra_t[j],&cm.metal_t[j],&cm.Rs_t[j],
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[3][j],&mk,&cm.cl_t[j],&cm.type_t[j]);
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or cm.metal_t[j]>0.06 || cm.cl_t[j]>5 || cm.type_t[j]>8.0 or 
    (cm.cl_t[j]==5 and float(cm.type_t[j])>8.0) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9) or (cm.cl_t[j]<5 and int(cm.type_t[j])==9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"mass: "<<cm.mass_t[j]<<"\t TefF:  "<<cm.Teff_t[j]<<"\t metal:  "<<cm.metal_t[j]<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;






////=================================== STELLAR HALO ============================================================ 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c.dat",'C','M','D','h');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&cm.age_h[j],&cm.logl_h[j],&cm.gra_h[j],&cm.metal_h[j],&cm.Rs_h[j],
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[3][j],&mk,&cm.cl_h[j],&cm.type_h[j]);
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || cm.age_h[j]<13.0 or cm.age_h[j]>15.0 or cm.cl_h[j]<0 or cm.cl_h[j]>5  or  cm.Teff_h[j]<0.0 or
    cm.metal_h[j]>0.04 || cm.cl_h[j]>8 || cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>7) or 
    (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
    cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
    cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;
  
}
///==============================================================//
///                                                              //
///                  optical_depth                               //
///                                                              //
///==============================================================//
void optical_depth(source & s)
{
    double ds =(double)s.nums*step;///kpc
    double CC=4.0*G*M_PI*ds*ds*pow(10.0,9.0)*M_sun/(velocity*velocity*KP);
    double dl,x,dx;
    s.od_disk=s.od_ThD=s.od_bulge=s.od_halo=s.opd=0.0;
    for(int k =1;k<s.nums;++k){
    dl =(double)k*step;///kpc
    x=dl/ds;
    dx=(double)step/ds/1.0;
    s.od_disk +=  s.rho_disk[k]*x*(1.0-x)*dx*CC;
    s.od_ThD +=   s.rho_ThD[k]*x*(1.0-x)*dx*CC;
    s.od_bulge += s.rho_bulge[k]*x*(1.0-x)*dx*CC;
    s.od_halo +=  s.rho_halo[k]*x*(1.0-x)*dx*CC;}
    s.opd= s.od_disk+s.od_ThD+s.od_bulge+s.od_halo;///total
    //cout<<"total_opticalD: "<<s.opd<<"\t od_disk: "<<s.od_disk<<endl;
   // cout<<"od_ThD: "<<s.od_ThD<<"\t od_bulge: "<<s.od_bulge<<"\t od_halo: "<<s.od_halo<<endl;
}
///==============================================================//
///                                                              //
///                  func_source   Initial amounts               //
///                                                              //
///==============================================================//
void func_source(source & s, CMD & cm, extinc & ex)
{
    int num,struc,nums;
    double rho,rf,Nblend[M]={0.0};
    double Akv,Avk,Alv;
    double Ds,Ai[M],Av;
    double Map[M];
   

    double maxnb=0.0;
    for(int i=0; i<M; ++i){
    s.Fluxb[i]=0.0;
    Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    Nblend[i]=Nblend[i]+RandN(sqrt(Nblend[i]),1.0);
    if(Nblend[i]<0.0)  Nblend[i]=0.0; 
    s.nsbl[i]=  double(Nblend[i]);
    if(Nblend[i]<=1.0) Nblend[i]=1.0;
    if(Nblend[i]>maxnb) maxnb=Nblend[i];}
    

    for(int k=1; k<=int(maxnb); ++k){
    do{
    num=int(fabs((double)rand()/(double)(RAND_MAX+1.)*Num*1.0));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[num] || num<5);///distance larger than 50.0 
    Ds=(double)(num*step);///in kpc
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}
   // cout<<"rho_disk:  "<<s.rho_disk[nums]<<"\t rho_bulge:  "<<s.rho_bulge[nums]<<"\t rho_ThD:  "<<s.rho_ThD[nums]<<"\t rho_halo:  "<<s.rho_halo[nums]<<endl;



    rf=fabs((double)rand()/(double)(RAND_MAX+1.))*fabs(s.Rostar0[nums] );
         if (rf<= fabs(s.rho_disk[nums]) ) struc=0;///thin disk
    else if (rf<=fabs(fabs(s.rho_disk[nums])+fabs(s.rho_bulge[nums]))) struc=1;/// bulge structure
    else if (rf<=fabs(fabs(s.rho_disk[nums])+fabs(s.rho_bulge[nums])+fabs(s.rho_ThD[nums]) ) ) struc=2;///thick disk
    else if (rf<=fabs(s.Rostar0[nums]) ) struc=3;///halo
    else {
    cout<<"Error rf:  "<<rf<<"\t  Rostar0[nums]:  "<<s.Rostar0[nums]<<"\t nums:  "<<nums<<endl;
    int uuee;  cin>>uuee; }
    
    if(k==1)    s.struc=struc;
   // cout<<"struc:  "<<struc<<endl;


   double Mab[M]; 
    if(struc==0){///thin disk
    
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){ 
    s.type=cm.type_d[num];
    s.mass=cm.mass_d[num];
    s.Teff=cm.Teff_d[num]; 
    s.col=Mab[2]-Mab[3]; 
    s.logl=cm.logl_d[num];
    s.cl=    cm.cl_d[num]; 
    s.logg= cm.gra_d[num]; 
    s.Rs=    cm.Rs_d[num];
    s.age=  cm.age_d[num]; }}


    if(struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    if(k==1) { 
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Teff=cm.Teff_b[num];
    s.col=Mab[2]-Mab[3]; 
    s.logl=cm.logl_b[num];
    s.cl=    cm.cl_b[num];
    s.logg =cm.gra_b[num]; 
    s.Rs =   cm.Rs_b[num];
    s.age=  cm.age_b[num]; }}


    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){ 
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Teff=cm.Teff_t[num];
    s.col=Mab[2]-Mab[3]; 
    s.logl=cm.logl_t[num];
    s.cl=    cm.cl_t[num];
    s.logg =cm.gra_t[num];
    s.Rs =   cm.Rs_t[num];
    s.age=  cm.age_t[num]; }}


    if(struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<M; ++i) { Mab[i]=cm.Mab_h[i][num];}
    if(k==1){
    s.cl  =  cm.cl_h[num]; 
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Teff=cm.Teff_h[num];
    s.col=Mab[2]-Mab[3]; 
    s.logl=cm.logl_h[num];
    s.logg =cm.gra_h[num]; 
    s.Rs =   cm.Rs_h[num];
    s.age=  cm.age_h[num];}}
    
    if(s.type>8.0 or s.type<2.0 or s.Rs<0.0 or s.Rs>1000.0 or s.mass<0.0 or 
    s.Teff<0.0 or s.cl==6 or s.logg<0.0 or s.logg<0.0  or s.age<0.0 or s.age>130.0 or Mab[0]<-10.0 or Mab[0]>30.0 or Mab[1]<-10.0 or Mab[1]>30.0 or Mab[2]<-10.0 or Mab[2]>30.0 or Mab[3]<-10.0 or Mab[3]>30.0){
     cout<<"ERROR:  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<endl; 
     cout<<"Mab_U:  "<<Mab[0]<<"\t Mab_B:  "<<Mab[1]<<"\t Mab_V:  "<<Mab[2]<<endl;
     cout<<"Mab_I:  "<<Mab[3]<<"\t num:  "<<num<<"\t struc:  "<<struc<<"\t Ds:  "<<Ds<<endl;
     int rre;  cin>>rre;}
    
    

   
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av>20.0 or Av<0.0 or Ds>20.0  or Ds<0.0){cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl; int yyw;  cin>>yyw; }
    if(Av<0.0)   Av=0.0;

    for(int i=0; i<M; ++i){    
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i],1.0);
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0) + Ai[i];
    if(Nblend[i]>=k)  s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));
    if(k==1){
    s.ext[i]=Ai[i];   s.Map[i]=Map[i];   s.Mab[i]= Mab[i];
    s.col=s.col+s.ext[2]-s.ext[3];}
    if(Ai[i]<0.0  or Ai[i]>100.0 or Map[i]<Mab[i] or Ds<0.0  or Ds>20.0  or Map[i]<0.0){
    cout<<"ERROR filter:  "<<i<<"\t extinction:  "<<Ai[i]<<"\t App_mag:  "<<Map[i]<<"\t Abso_mag:  "<<Mab[i]<<endl;
    int rre; cin>>rre;}}
    }///loop 


    for(int i=0; i<M; ++i){
    if(s.Fluxb[i]<=0.0){cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<endl; int yye; cin>>yye; }
    s.magb[i] = -2.5*log10( s.Fluxb[i] );
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    if(int(s.nsbl[i])<0 or (Nblend[i]==1.0 && s.blend[i]<1.0) or s.blend[i]>1.0  or s.blend[i]<0.0){ 
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t Nlend: "<<Nblend[i]<<"\t s.blend  "<<s.blend[i]<<endl; int uue; cin>>uue;}}
 
 
   // cout<<"Ds:  "<<s.Ds<<"\t nums:  "<<s.nums<<endl;
   // cout<<"End of func source "<<endl;
   // s.magni=fabs((s.u0*s.u0+2.0)/(s.u0*sqrt(s.u0*s.u0+4.0)));
  // cout<<"Map[0]:  "<<s.Map[0]<<"\t Mab[0]:  "<<s.Mab[0]<<"\t extinc:  "<<s.ext[0]<<endl;
  // cout<<"Map[3]:  "<<s.Map[3]<<"\t Mab[3]:  "<<s.Mab[3]<<"\t extinc:  "<<s.ext[3]<<endl;
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
   //  cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }
    // if(ex.dis[i]==8.25) ex.exks8= ex.Extks[i]; 
     }}
     cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
///==============================================================//
///                                                              //
///                  func_lens     Initial amounts               //
///                                                              //
///==============================================================//
void func_lens(lens & l, source & s, CMD & cm, extinc & ex)
{
    double f,test;
    double rholens[s.nums+2]={0.0};
    l.rhomaxl=0.0;

    double f1, f2, hel, minm, dism; 
    int il; 
    double Akv,Avk,Alv;
    double Ai[M],Av;
    double Mabl[M]; 


    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds) {cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}


    do{
    l.numl = (int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test = ((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;
    int ue; cin>>ue;}
    }while(test>rholens[l.numl]);


   double  randflag= ((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
       if (randflag<= s.rho_disk[l.numl]) l.struc=0;///thin disk
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
  else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
  else if (randflag<=  s.Rostar0[l.numl]) l.struc=3;//halo
  else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}

///#################################################################################
  double mp= double(13.0*0.000954588); 
  double mmin=double(1.0*0.000954588); 
  if(l.struc==0){///thin disk
  
  
  
  f1=pow(mmin,-1.0);
  f2=pow(1.0,-1.6)*pow(0.08,-0.7+1.6)*pow(mp,-1.0+0.7); 
  do{
  l.Ml=((double)rand()/(double)(RAND_MAX+1.)  )*(1.0-mmin) + mmin;
  test=((double)rand()/(double)(RAND_MAX+1.)  )*(f1-f2) +f2 ;
  if(l.Ml<mp)           f=pow(l.Ml,-1.0);///planet
  else if(l.Ml<=0.08)   f=pow(l.Ml,-0.7)*pow(mp,-1.0+0.7);//BD 
  else if(l.Ml<=1.0 )   f=pow(l.Ml,-1.6)*pow(0.08,-0.7+1.6)*pow(mp,-1.0+0.7);
  }while(test>f);
  
  
  if(l.Ml>0.08){
  minm=10000.0, dism=0.0; il=-1;
  for(int i=0; i<N1;  ++i){
  dism=fabs(cm.mass_d[i]-l.Ml);
  if(dism<minm and cm.cl_d[i]>2 and cm.cl_d[i]<7){il=i;    minm=dism;}}
  if(il==-1 or minm==10000 or cm.cl_d[il]<=2 or cm.cl_d[il]>=7){
  cout<<"Error il: "<<il<<"\t minm: "<<minm<<"\t cl: "<<cm.cl_d[il]<<endl;  int ooi;  cin>>ooi; }
  for(int i=0; i<M; ++i){Mabl[i]=cm.Mab_d[i][il];}
  l.cl =   cm.cl_d[il];  
  l.type=cm.type_d[il];
  l.Teff=cm.Teff_d[il]; 
  l.logl=cm.logl_d[il];
  l.logg= cm.gra_d[il]; 
  l.Rl=    cm.Rs_d[il]; 
  l.age=  cm.age_d[il]; 
  if(l.Ml<0.0 or l.cl<0 or l.type<00.0 or l.Teff<0.0 or l.Rl<0.0 or l.logg<0.0 or l.age<0.0){
  cout<<"Error il: "<<il<<"\t mass: "<<l.Ml<<endl;  int yye; cin>>yye;}}
  else if(l.Ml>mp) l.type=9;
  else             l.type=10.0; }
  
  
  
///#################################################################################
  if(l.struc==1){///Galactic bulge
  
  f1=pow(mmin,-1.0);
  f2=pow(1.0,-2.35)*pow(0.08,-0.7+2.35)*pow(mp,-1.0+0.7); 
  do{
  l.Ml=((double)rand()/(double)(RAND_MAX+1.)  )*(1.0-mmin) + mmin;
  test=((double)rand()/(double)(RAND_MAX+1.)  )*(f1-f2) +f2 ;
  if(l.Ml<mp)           f=pow(l.Ml,-1.0);///planet
  else if(l.Ml<=0.08)   f=pow(l.Ml,-0.7)*pow(mp, -1.0+0.7);//BD 
  else if(l.Ml>=0.08)   f=pow(l.Ml,-2.35)*pow(0.08,-0.7+2.35)*pow(mp,-1.0+0.7); //
  }while(test>f);
    
    
  if(l.Ml>0.08){
  minm=10000.0, dism=0.0; il=-1;
  for(int i=0; i<N2;  ++i){
  dism=fabs(cm.mass_b[i]-l.Ml);
  if(dism<minm and cm.cl_b[i]>2 and cm.cl_b[i]<7){ il=i;   minm=dism; }}
  if(il==-1 or minm==10000 or cm.cl_b[il]<=2 or cm.cl_b[il]>=7){
  cout<<"Error il: "<<il<<"\t minm: "<<minm<<"\t cl: "<<cm.cl_b[il]<<endl;  int ooi;  cin>>ooi; }
  for(int i=0; i<M; ++i){Mabl[i]=cm.Mab_b[i][il];}
  l.cl  =  cm.cl_b[il];  
  l.type=cm.type_b[il];
  l.Teff=cm.Teff_b[il]; 
  l.logl=cm.logl_b[il];
  l.logg= cm.gra_b[il]; 
  l.Rl=    cm.Rs_b[il]; 
  l.age=  cm.age_b[il]; 
  if(l.Ml<0.0 or cm.cl_b[il]<0 or l.type<00.0 or l.Teff<0.0 or l.Rl<0.0 or l.age<0.0){
  cout<<"Error mass_b: "<<cm.mass_b[il]<<"\t mass: "<<s.mass<<endl;  int yye; cin>>yye; }}
  else if(l.Ml>mp) l.type=9;
  else             l.type=10.0; }



///#################################################################################
  if(l.struc==2){///thick disk 
  
  f1=pow(mmin,-1.0);
  f2=pow(1.0,-0.5)*pow(0.08,-0.7+0.5)*pow(mp,-1.0+0.7); 
  do{
  l.Ml=((double)rand()/(double)(RAND_MAX+1.)  )*(1.0-mmin) + mmin;
  test=((double)rand()/(double)(RAND_MAX+1.)  )*(f1-f2) +f2 ;
  if(l.Ml<mp)           f=pow(l.Ml,-1.0);///planet
  else if(l.Ml<=0.08)   f=pow(l.Ml,-0.7)*pow(mp, -1.0+0.7);//BD 
  else if(l.Ml>=0.08)   f=pow(l.Ml,-0.5)*pow(0.08,-0.7+0.5)*pow(mp,-1.0+0.7); //
  }while(test>f);
   
  if(l.Ml>0.08){
  minm=10000.0, dism=0.0; il=-1;
  for(int i=0; i<N3;  ++i){
  dism=fabs(cm.mass_t[i]-l.Ml);
  if(dism<minm and cm.cl_t[i]>2 and cm.cl_t[i]<7){ il=i;   minm=dism; }}
  if(il==-1 or minm==10000 or cm.cl_t[il]<=2 or cm.cl_t[il]>=7){
  cout<<"Error il: "<<il<<"\t minm: "<<minm<<"\t cl: "<<cm.cl_t[il]<<endl;  int ooi;  cin>>ooi; }
  
  for(int i=0; i<M; ++i){Mabl[i]=cm.Mab_t[i][il];}
  l.cl=    cm.cl_t[il];  
  l.type=cm.type_t[il];
  l.Teff=cm.Teff_t[il]; 
  l.logl=cm.logl_t[il];
  l.logg= cm.gra_t[il]; 
  l.Rl=    cm.Rs_t[il]; 
  l.age=  cm.age_t[il]; 
  if(l.Ml<0.0 or cm.cl_t[il]<0 or l.type<00.0 or l.Teff<0.0 or l.Rl<0.0 or l.age<0.0){
  cout<<"Error mass_t: "<<cm.mass_t[il]<<"\t mass: "<<s.mass<<endl;  int yye; cin>>yye; } }
  else if(l.Ml>mp) l.type=9;
  else             l.type=10.0;}

 
 
 
///#################################################################################
  if(l.struc==3){///stellar halo
  
  f1=pow(mmin,-1.0);
  f2=pow(1.0,-0.5)*pow(0.08,-0.7+0.5)*pow(mp,-1.0+0.7);  
  do{
  l.Ml=((double)rand()/(double)(RAND_MAX+1.)  )*(1.0-mmin) + mmin;
  test=((double)rand()/(double)(RAND_MAX+1.)  )*(f1-f2) +f2 ;
  if(l.Ml<mp)           f=pow(l.Ml,-1.0);///planet
  else if(l.Ml<=0.08)   f=pow(l.Ml,-0.7)*pow(mp, -1.0+0.7);//BD 
  else if(l.Ml>=0.08)   f=pow(l.Ml,-0.5)*pow(0.08,-0.7+0.5)*pow(mp,-1.0+0.7); //
  }while(test>f);
  
  if(l.Ml>0.08){
  double minm=10000.0, dism=0.0; int il=-1;
  for(int i=0; i<N3;  ++i){
  dism=fabs(cm.mass_h[i]-l.Ml);
  if(dism<minm and cm.cl_h[i]>2 and cm.cl_h[i]<7){ il=i;   minm=dism; }}
  if(il==-1 or minm==10000 or cm.cl_h[il]<=2 or cm.cl_h[il]>=7) {
  cout<<"Error il: "<<il<<"\t minm: "<<minm<<"\t cl: "<<cm.cl_h[il]<<endl;  int ooi;  cin>>ooi; }
  
  for(int i=0; i<M; ++i){Mabl[i]=cm.Mab_h[i][il];}
  l.cl=    cm.cl_h[il];  
  l.type=cm.type_h[il];
  l.Teff=cm.Teff_h[il]; 
  l.logl=cm.logl_h[il];
  l.logg= cm.gra_h[il]; 
  l.Rl=    cm.Rs_h[il]; 
  l.age=  cm.age_h[il]; 
  if(l.Ml<0.0 or cm.cl_h[il]<0 or l.type<00.0 or l.Teff<0.0 or l.Rl<0.0 or l.age<0.0){
  cout<<"Error mass_h: "<<cm.mass_h[il]<<"\t mass: "<<s.mass<<endl;  int yye; cin>>yye;} }
  else if(l.Ml>mp) l.type=9;
  else             l.type=10.0;}


  l.Dl=l.numl*step;///kpc
  l.xls=l.Dl/s.Ds;
  l.RE=sqrt(4.0*G*l.Ml*M_sun*s.Ds*KP)/velocity;
  l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
  vrel(s,l);
  l.tE=l.RE/(l.Vt*1000.0*3600.0*24.0); 
  s.ro_star=double(s.Rs*Rsun*l.xls/l.RE); 
  l.mul=  l.Vt*1000.0*180.0*3600.0*1000.0*3600.0*24.0*365.24222/(l.Dl*KP*M_PI);
  l.u0=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0003745972354);
  l.pt1=-1.5*l.tE;  
  l.pt2= 1.5*l.tE;  
  l.dt=double(10.0/60.0/24.0);  ///10 min
 
}
///===========================================================================//
double RandN(double sigma, double nnd){
   double rrd,f,frand;
   do{
   rrd=double(((double)rand()/(double)(RAND_MAX+1.))*2.0-1.0)*sigma*nnd; ///[-N sigma:N sigma]
   f= exp(-0.5*rrd*rrd/sigma/sigma);
   frand=fabs((double)rand()/((double)(RAND_MAX+1.))*1.0);
   }while(frand>f);
   return(rrd);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars  
   double fh=1.0;///No limitation 
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;

   for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = R_sun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/R_sun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/R_sun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[M_sun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[M_sun/deg^2]
s.Nstari[i]= binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3] 

s.Nstari[i]=s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]); 
}
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
 if (l.Dl==0.0) l.Dl=0.00034735;
  double pi=M_PI;
  double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
  double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
  if(Rlc==0.0) Rlc=0.0000000000034346123;
  if(Rsc==0.0) Rsc=0.000000000004762654134; 
 
  double LVx, SVx;
  double SVT, SVR, SVZ, LVT, LVR, LVZ;
  double fv, testfv, test, age;
  double  VSunx, vls2, vls1;
  double betal, betas, deltal, deltas, tetd ;



  double NN=2.5;
  double sigma_R_Disk,   sigma_T_Disk,  sigma_Z_Disk;
  double sigma_R_DiskL,  sigma_T_DiskL, sigma_Z_DiskL;
  double sigma_R_DiskS,  sigma_T_DiskS, sigma_Z_DiskS;
  double sigma_R_TDisk=67.0,  sigma_T_TDisk=51.0, sigma_Z_TDisk=42.0;
  double sigma_R_halo= 131.0, sigma_T_halo=106.0, sigma_Z_halo=85.0;
  double sigma_R_Bulge=113.0, sigma_T_Bulge=115.0, sigma_Z_Bulge=100.0;
  double Rho[8]={00.0}; double maxr=0.0;
  for(int i=0; i<8; ++i){Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}

 
for (int i=0;i<2; ++i){
 test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])                       {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0; age= 0.075;}
else if(test<=(Rho[0]+Rho[1]))              {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0; age=0.575; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]))       {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;age=1.5;  }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3])){sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2; age=2.5; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))       {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8; age=4.0; }
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5])){sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4; age=6.0; }
else if(test<=maxr)                                       {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5; age=8.5; }
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}
    if(i==0) {
    sigma_R_DiskS= sigma_R_Disk;
    sigma_T_DiskS= sigma_T_Disk;
    sigma_Z_DiskS= sigma_Z_Disk;}
    if(i==1){
    sigma_R_DiskL= sigma_R_Disk;
    sigma_T_DiskL= sigma_T_Disk;
    sigma_Z_DiskL= sigma_Z_Disk; }}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.struc==0){///Galactic disk
    SVR= RandN(sigma_R_DiskS, NN);
    SVT= RandN(sigma_T_DiskS, NN);
    SVZ= RandN(sigma_Z_DiskS, NN); }

    else if(s.struc==1){///Galactic bulge
    SVR= RandN(sigma_R_Bulge, NN);
    SVT= RandN(sigma_T_Bulge, NN);
    SVZ= RandN(sigma_Z_Bulge, NN); }

    else if(s.struc==2){///thick disk
    SVR= RandN(sigma_R_TDisk, NN);
    SVT= RandN(sigma_T_TDisk, NN);
    SVZ= RandN(sigma_Z_TDisk, NN); }

    else if(s.struc==3){///stellar halo
    SVR= RandN(sigma_R_halo, NN);
    SVT= RandN(sigma_T_halo, NN);
    SVZ= RandN(sigma_Z_halo, NN); }
    if(s.struc==0 or s.struc==2)  SVT =SVT+ vro_sun*(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712);
    l.vs=sqrt( SVR*SVR + SVT*SVT + SVZ*SVZ );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
   if(l.struc==0){///Galactic disk
   LVR= RandN(sigma_R_DiskL, NN) ;
   LVT= RandN(sigma_T_DiskL, NN) ;
   LVZ= RandN(sigma_Z_DiskL, NN) ; }

   else if(l.struc==1){///Galactic bulge
   LVR= RandN(sigma_R_Bulge, NN) ;
   LVT= RandN(sigma_T_Bulge, NN) ;
   LVZ= RandN(sigma_Z_Bulge, NN) ; }

   else if(l.struc==2){///thick disk
   LVR= RandN(sigma_R_TDisk, NN) ;
   LVT= RandN(sigma_T_TDisk, NN) ;
   LVZ= RandN(sigma_Z_TDisk, NN) ; }

   else if(l.struc==3){///stellar halo
   LVR= RandN(sigma_R_halo, NN);
   LVT= RandN(sigma_T_halo, NN);
   LVZ= RandN(sigma_Z_halo, NN); }
   if(l.struc==0 or l.struc==2)  LVT = LVT+vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712);
   l.vl=sqrt( LVT*LVT + LVZ*LVZ + LVR*LVR );
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH  BETA  HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    betal=betas=0.0;
    tetd= s.TET;
    test= double(l.Dl*cos(s.FI)*sin(tetd)/Rlc);
    if(fabs(test-1.0)<0.01)       betal= pi/2.0;
    else if(fabs(test+1.0)<0.01)  betal=-pi/2.0;
    else                          betal=asin(test);
    
    test= double(s.Ds*cos(s.FI)*sin(tetd)/Rsc); 
    if( fabs(test-1.0)<0.01)     betas=pi/2.0;
    else if(fabs(test+1.0)<0.01) betas=-pi/2.0;
    else                         betas=asin(test);
    
    if(R_sun < fabs(l.Dl*cos(s.FI)*cos(tetd))) betal= pi-betal; 
    if(R_sun < fabs(s.Ds*cos(s.FI)*cos(tetd))) betas= pi-betas; 

   if(fabs(l.Dl*cos(s.FI)*sin(tetd))>Rlc or fabs(test)>1.0){
   cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
   cout<<"FI: "<<s.FI<<"\t TET: "<<tetd<<"\t betal:  "<<betal<<endl;
   cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<"\t betas:   "<<betas<<endl;
   cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(tetd)/Rlc<<"\t sin(s): "<<test<<endl;
   int ew; cin>>ew;}
///HHHHHHHHHHHHHHHHHHHHHHHHHH  DELTA   HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    if(s.TET>pi)  tetd=s.TET-2.0*pi; 
    deltal= pi - fabs(tetd) -fabs(betal);
    deltas= pi - fabs(tetd) -fabs(betas);  
    if(betal<0.0)  deltal= -1.0*deltal;
    if(betas<0.0)  deltas= -1.0*deltas;
    s.deltao= double( pi-fabs(tetd) );
    if(tetd<0.0)  s.deltao=-1.0*s.deltao; 

///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    s.SV_n1 =+SVR * sin(deltas)- SVT * cos(deltas);
    s.LV_n1 =+LVR * sin(deltal)- LVT * cos(deltal);
    s.VSun_n1=+VSunR*sin(s.deltao)-VSunT*cos(s.deltao);
    
    SVx= -SVR*cos(deltas)- SVT*sin(deltas);
    LVx= -LVR*cos(deltal)- LVT*sin(deltal);
    VSunx= -VSunR*cos(s.deltao) -VSunT*sin(s.deltao);
    
    s.SV_n2=-sin(s.FI)*(SVx) + cos(s.FI)*SVZ;
    s.LV_n2=-sin(s.FI)*(LVx) + cos(s.FI)*LVZ;
    s.VSun_n2=-sin(s.FI)*(VSunx)+cos(s.FI)*(VSunZ);
 
    
    vls1= l.xls*s.SV_n1 - s.LV_n1 +(1.0-l.xls)*s.VSun_n1;  ///Source - lens 
    vls2= l.xls*s.SV_n2 - s.LV_n2 +(1.0-l.xls)*s.VSun_n2;  /// Source -lens
    l.Vt=sqrt(fabs( vls1*vls1 + vls2*vls2 ) );
    
    
    if(vls1==0.0 and vls2==0.0){
    cout<<"Big error both zeros:  "<<vls1<<"\t vls_b:  "<<vls2<<endl;
    int iee;  cin>>iee;}

    if (l.Vt<0.0 or l.Vt>1.0e6){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
    cout<<"source type:  "<<s.type<<"\t s.struc:  "<<s.struc<<endl;
    int yee; cin>>yee;}
    //cout<<"(vrel_func):   s.deltao:  "<<s.deltao<<"FI: "<<s.FI<<endl;
}
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

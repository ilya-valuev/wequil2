#include "analyse.h"
#include "interact.h"
                                                     
#include "four.h"
#include "common.h"

# ifdef  UNIX

#include <sys/stat.h> // define mkdir only
#include <sys/types.h>
#include <unistd.h> // defines execl only

#include <errno.h> // defines errno only

# endif


/*
# ifndef log2
double log2(double x){ return log(x)/log(2); }
# endif */

int xx_used=0;
int time_aver=1;

char out_dir[200]="";

double random2(void){
 return 1.-2.*((double)rand())/RAND_MAX;
}

double anL,anGamma,anT,anRDebye;
double animass;
int npart,nion,nelec;
Vector_3 *anx, *anv;
Vector_3 flow[2];
double dt;
long nsteps,cur_step;
double ttime;
double Wp;


double amod12(double a, double b){
  return a-((long)(a/b)*b);
}


// returns vector between -L/2 and L/2
Vector_3 r_cell(Vector_3 r){
 int i;
 Vector_3 r1;
 for(i=0;i<3;i++){
   if(r[i]>0)r1[i]=amod12(r[i]+anL/2,anL)-anL/2;
   else r1[i]=amod12(r[i]-anL/2,anL)+anL/2;
 }
 return r1;
}



long gettypesize(long t, int n){
  long sz=0;
  if(t== NOTHING)return 0;
  if(t&COORD)sz+=3*n*sizeof(double);
  if(t&VEL)sz+=3*n*sizeof(double);
  if(t&FLOW)sz+=2*3*sizeof(double);
  return sz;
}



 


FILE *tFile;
long tType=0;
long rnsteps;

int InitTrajectory(long &type,long step, char *datafile,char *dname,double &l,double &g,double &t,double &imass){
  tFile=Err_fopen(datafile,"rb");
  char buffer[200];
  fread(buffer,sizeof(char),200,tFile); // skipping text header
  fread(buffer,sizeof(char),200,tFile); // reading
  int nn, res;  
  imass=1.;
  res=sscanf(buffer,"%d %lf %lf %lf %s %lf %ld %lf %ld %lf", &nn, &l,&g,&t,dname,&dt,&nsteps,&Wp,&type,&imass);
  if(res==8)type=VEL|COORD;
  else if(res!=9 && res !=10)serror("Invalid trajectory file: %s\n",datafile);


  tType=type;
  fseek(tFile,0,SEEK_END);
  double lp=(double)ftell(tFile);
  lp-=400;
  long full=gettypesize(COORD|VEL|type,nn);
  lp-=full;

  long dsz=gettypesize(type,nn);
  if(!dsz)lp=0;
  else lp/=dsz;
  //lp/=4;
  //lp/=6;
  //lp/=nn;
  if(((long)lp)!=nsteps-1){
    printf("Warning: incomplete trajectory file\n"
           "reported steps: %ld,  actual steps: %.1f\n",nsteps,lp+1);
    nsteps=(long)lp+1;
    printf("Setting nsteps to %ld\n",nsteps);
  }
  printf("Trajectory length: %f\n",nsteps*dt*Wp);

  if(step>=0)fseek(tFile,400+dsz*step,SEEK_SET);
  else{
   fseek(tFile,400+dsz*(nsteps-1),SEEK_SET);
   dsz=gettypesize(type&(~(VEL|COORD)),nn);
   fseek(tFile,dsz,SEEK_CUR);
  }

  ttime=0.;
  cur_step=0;
  rnsteps=nsteps;
  return nn;
}
  

void EndTrajectory(){
 fclose(tFile);
 ttime=0;
 nsteps=cur_step=0;
}

extern long ref_type;
extern long rnsteps;
extern int trjstep;
int BadStatus=0;

double ReadStep(long type, int n){
  static long count =0;
  BadStatus=1;
  int i;
  xx_used=0;
  for(i=0;i<trjstep;i++){
    count++;
    if(count==rnsteps-1){
      type|=(COORD|VEL);
      ref_type=tType=type;
    }



    if(count>=rnsteps){
      //cur_step--;
      fseek(tFile,400,SEEK_SET);
      return -1.;
    }
    if(feof(tFile))return -2.;
    if(tType&type&FLOW)fread(flow,gettypesize(FLOW,n),1,tFile);
    else if(tType&FLOW)fseek(tFile,gettypesize(FLOW,n),SEEK_CUR);

    if(feof(tFile))return -2.;
    if(tType&type&COORD)fread(anx,sizeof(Vector_3),n,tFile);
    else if(tType&COORD)fseek(tFile,gettypesize(COORD,n),SEEK_CUR);

    if(feof(tFile))return -2.;
    if(tType&type&VEL)fread(anv,sizeof(Vector_3),n,tFile);
    else if(tType&VEL)fseek(tFile,gettypesize(COORD,n),SEEK_CUR);
  }

  cur_step++;
  ttime+=dt;
  BadStatus=0;
  return ttime;
}


Correlation *vc_Ci, *vc_Ce;



int vc_ion=0, vc_elc=0;
int vc_type=VC_FFT;
int vc_n, vc_sav;

void Vcorr_init(){
  if(!(tType&VEL))serror("Can't calculate velocity autocorrelation,\n"
			   " trajectory file does not contain Velocity data.\n");
  Set_stop(1);
  int i;
  char tmpstr[256];

  Read_param("Ion velocity acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_ion=1;
  else vc_ion=0;

  Read_param("Electron velocity acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_elc=1;
  else vc_elc=0;


  Read_param("Substract mean: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_sav=1;
  else vc_sav=0;

  Read_param("Calculation type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))vc_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))vc_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))vc_type=VC_STRAIT;
  else serror("Velocity acf: unknown calculation type.\n");

  if(vc_type==VC_FFT){
    printf("Velocity acf: fourier degree: %d",(int)log2(nsteps)+1);
    vc_n=1<<((int)log2(nsteps)+1);
    printf(", %d points\n",vc_n);
  }
  else vc_n=nsteps;

  if(vc_ion){
   vc_Ci=new Correlation[3*nion];
   if(!vc_Ci)serror("Velocity acf: memory allocation error.\n");
   for(i=0;i<3*nion;i++){
    vc_Ci[i].init(vc_n,nsteps*dt);
    vc_Ci[i].sub_mean=vc_sav;
    vc_Ci[i].insert_begin(nsteps);
   }
  }

  if(vc_elc){
   vc_Ce=new Correlation[3*nelec];
   if(!vc_Ce)serror("Velocity acf: memory allocation error.\n");
   for(i=0;i<3*nelec;i++){
    vc_Ce[i].init(vc_n,nsteps*dt);
    vc_Ce[i].sub_mean=vc_sav;
    vc_Ce[i].insert_begin(nsteps);
   }
  }
}


void Vcorr_step(){
  int i,j;
  if(vc_ion){
   for(i=0;i<nion;i++){
    for(j=0;j<3;j++){
     vc_Ci[3*i+j].insert_next(anv[i][j]);
    }
   }
  }
  if(vc_elc){
   for(i=0;i<nelec;i++){
    for(j=0;j<3;j++){
     vc_Ce[3*i+j].insert_next(anv[nion+i][j]);
    }
   }
  }

}

void output_vcorr(char *filename,Correlation* Corr, int npart){
  int i,j;

  double *acf_data;

  if(vc_type==VC_DIRECT){
   strcat(filename,"d.dat");
   printf("Direct correlation ...\n");
   for(i=0;i<npart;i++){
    //printf("Direct correlation #%d...\n",i+1);
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
     Corr[3*i+j].direct();
    }
   }
   for(j=0;j<Corr[0].n;j++){
    //Corr[0].ini[j]=Corr[0].arr[j]*Corr[0].arr[j];
    for(i=1;i<3*npart;i++){
     Corr[0].arr[j]+=Corr[i].arr[j];
     //Corr[0].ini[j]+=Corr[i].arr[j]*Corr[i].arr[j]; // average square
    }
    Corr[0].arr[j]/=(double)(3*npart);
    //Corr[0].ini[j]/=(double)(3*npart);
    //Corr[0].ini[j]=Corr[0].ini[j]-Corr[0].arr[j]*Corr[0].arr[j];
   }
   Corr[0].normalize();
   //x*=x;
   //for(j=0;j<Corr[0].n;j++)Corr[0].ini[j]/=x;
   acf_data=Corr[0].arr;
  }
  else if(vc_type==VC_STRAIT){
   strcat(filename,"s.dat");
   printf("Strait correlation ...\n");
   for(i=0;i<npart;i++){
    //printf("Strait correlation #%d...\n",i+1);
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
     Corr[3*i+j].strait();
    }
   }
   for(j=0;j<Corr[0].n;j++){
    //Corr[0].ini[j]=Corr[0].arr[j]*Corr[0].arr[j];
    for(i=1;i<3*npart;i++){
     Corr[0].arr[j]+=Corr[i].arr[j];
     //Corr[0].ini[j]+=Corr[i].arr[j]*Corr[i].arr[j]; // average square
    }
    Corr[0].arr[j]/=(double)(3*npart);
    //Corr[0].ini[j]/=(double)(3*npart);
    //Corr[0].ini[j]=Corr[0].ini[j]-Corr[0].arr[j]*Corr[0].arr[j];
   }
   Corr[0].normalize();
   //x*=x;
   //for(j=0;j<Corr[0].n;j++)Corr[0].ini[j]/=x;
   acf_data=Corr[0].arr;
  }
  else{
   strcat(filename,"f.dat");

   for(i=0;i<npart;i++){
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
    }
   }

   printf("Velocity auto corelations: performing  FFT...\n");
   acf_data=Calculate_forth(Corr[0].n,Corr,3*npart);
   printf("Inverse FFT...\n");
   Calculate_back(acf_data,Corr[0].n);
   double c=acf_data[0];
   for(i=0;i<Corr[0].n;i++)acf_data[i]=acf_data[2*i]/c;
  }
  TableFunction acf(Corr[0].n,acf_data);

  double tmax=Corr[0].maxtime*Wp;
  acf.xscale(0.,tmax);
  //Corr[0].fini.xscale(0.,tmax);
  FILE *f1=Err_fopen(filename,"wt");
  //fprintf(f1,"#1-t 2-vcorr 3-d^2\n");
  double dt=0.05,t;
  for(t=0.;t<tmax /*&& t< 150.*/;t+=dt){
    //if(t>8.15){
    // printf("!");
    //}
    fprintf(f1,"%f  %f",t,acf(t));
    //if(vc_type==VC_STRAIT || vc_type==VC_DIRECT){
    // fprintf(f1," %f",Corr[0].fini(t));
    //}
    fprintf(f1,"\n");
  }
  fclose(f1);
}

void Vcorr_process(char *dname){
  char str[256]="vcorr.dat";

  if(vc_ion){
   strcpy(str,out_dir);
   strcat(str,dname);
   strcat(str,"-vi");
   output_vcorr(str,vc_Ci,nion);
   delete[] vc_Ci;
  }
  if(vc_elc){
   strcpy(str,out_dir);
   strcat(str,dname);
   strcat(str,"-ve");
   output_vcorr(str,vc_Ce,nelec);
   delete[] vc_Ce;
  }
}



// DM calculation
int dm_sav=0, dm_type, dm_n;
Correlation dm_C[3];
//Correlation *dm_C;

void DM_init(){
  if(!(tType&COORD))serror("Can't calculate dipole moment,\n"
			   " trajectory file does not contain Coords data.\n");
  Set_stop(1);
  int i;
  char tmpstr[256];

  Read_param("DM substract mean: %s",tmpstr);
  if(strstr(tmpstr,"y"))dm_sav=1;
  else dm_sav=0;

  Read_param("DM acf type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))dm_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))dm_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))dm_type=VC_STRAIT;
  else serror("DM acf: unknown calculation type.\n");

  if(dm_type==VC_FFT){
    printf("DM acf: fourier degree: %d",(int)log2(nsteps)+1);
    dm_n=1<<((int)log2(nsteps)+1);
    printf(", %d points\n",dm_n);
  }
  else dm_n=nsteps;

  for(i=0;i<3;i++){
   dm_C[i].init(dm_n,nsteps*dt*Wp);
   dm_C[i].sub_mean=dm_sav;
   dm_C[i].insert_begin(nsteps);
  }
}


void DM_step(){
  int i;
  Vector_3 moment(0.,0.,0.);
  //printf("istart\n");

  for(i=0;i<nion;i++){
   if(!xx_used)moment+=r_cell(anx[i]);
   else moment+=anx[i];
  }

  for(i=0;i<nelec;i++){
   if(!xx_used)moment-=r_cell(anx[nion+i]);
   else moment-=anx[nion];
  }
  for(i=0;i<3;i++){
   dm_C[i].insert_next(moment[i]);
  }
  //printf("iend\n");
}


void DM_process(char *dname){

  char filename[256]="dm.dat";


  strcpy(filename,out_dir);
  strcat(filename,dname);
  strcat(filename,"-dm");

  int i,j;


  double *acf_data;

  for(i=0;i<3;i++)dm_C[i].insert_end();

  if(dm_type==VC_DIRECT || dm_type==VC_STRAIT){
   if(dm_type==VC_DIRECT){
    //strcat(filename,"d.dat");
    strcat(filename,".dat");
    printf("Dipole moment: Direct correlation ...\n");
    for(i=0;i<3;i++){
      dm_C[i].direct();
    }
   }
   else{
    //strcat(filename,"s.dat");
    strcat(filename,".dat");
    printf("Dipole moment: Strait correlation ...\n");
    for(j=0;j<3;j++){
      dm_C[j].strait();
    }
   }
   for(j=0;j<dm_C[0].n;j++){
    for(i=1;i<3;i++){
     dm_C[0].arr[j]+=dm_C[i].arr[j];
     dm_C[0].arr[j+dm_C[0].n]+=dm_C[i].arr[j+dm_C[0].n];
    }
    dm_C[0].arr[j]/=3;
    dm_C[0].arr[j+dm_C[0].n]/=3;
   }

  }
  else{ // FFT
   //strcat(filename,"f.dat");
   strcat(filename,".dat");

   printf("Dipole moment autocorelation: performing  FFT...\n");

   acf_data=Calculate_forth(dm_C[0].n,dm_C,3);
   printf("Inverse FFT...\n");
   Calculate_back(acf_data,dm_C[0].n);

   for(i=0;i<dm_C[0].n;i++){
    dm_C[0].arr[i]=acf_data[2*i];
    dm_C[0].arr[i+dm_C[0].n]=acf_data[2*dm_C[0].n+i];
   }
   delete acf_data;
  }
  dm_C[0].normalize();

  double tmax=dm_C[0].maxtime;

  for(i=0;i<3;i++)dm_C[i].fini.xscale(0.,tmax);

  FILE *f1=Err_fopen(filename,"wt");
  fprintf(f1,"#1-t 2-Dx 3-Dy 4-Dz 5-D^2(t) 6-acf(D) 7-w 8-Fourier^2(acf(D))\n");

  double dt=0.1,t,w=0.,dw=(2*M_PI*dm_C[0].n/tmax)*(dt/tmax);
  Vector_3 d;
  for(t=0.;t<tmax /*&& t< 150.*/;t+=dt){
    for(i=0;i<3;i++)d[i]=dm_C[i].data(t);

    fprintf(f1,"%f %f %f %f %f %f",t,d[0],d[1],d[2],d*d,dm_C[0].r(t));
    if(dm_type==VC_FFT){
     fprintf(f1, " %f %f\n",w,dm_C[0].four_sq(w));
    }
    else fprintf(f1,"\n");
    w+=dw;
  }
  fclose(f1);
}







int mc_np=0, mc_rnd=0, mc_nop=0;
Vector_3 *mc_tpx;
char **mc_tpname;
double mc_r0, mc_r1;
int mc_nr, mc_xmgr=0, mc_prcoord=0, mc_prav,mc_correl=0;
cList* mc_plist=NULL;
double mc_F0;
int mc_type, mc_force=0x1, mc_ntype=0 /* ion*/;

enum {ON_SPACE,ON_ION,ON_ELC};
enum {COULOMB=0x1, EFFECTIVE=0x2, ALL=0x3};

typedef Distribution Distrvect[4];
//typedef Correlation  Corrvect[3];
Distrvect *Mf;
Distrvect *Ef;
Correlation *MfCorr;
Correlation *EfCorr;
//Distrvect Mfav;

int micro_field=0;

EffPotential mc_potential;


void Microfield_init(){
  if(!(tType&COORD))serror("Can't calculate microfield,\n"
			   " trajectory file does not contain Coords data.\n");
  Set_stop(1);
  int i,j;
  char tmpstr[256];
  int ac_n=1;

  Read_param("Microfield at: %s",tmpstr);
  if(strstr(tmpstr,"space")){
    micro_field=2;
    int &np=mc_np;
    Read_param("Average points: %d %s", &np, tmpstr);
    mc_type=ON_SPACE;
    mc_tpx= new Vector_3[np];
  }
  else if(strstr(tmpstr,"ion")){
    mc_np=nion;
    mc_type=ON_ION;
    mc_tpx=anx+nelec;
    
  }
  else if(strstr(tmpstr,"electron")){
    mc_np=nelec;
    mc_type=ON_ELC;
    mc_tpx=anx;
    
  }
    
  else serror("Unrecognized microfield calculation type: %s\n",tmpstr);
 
  typedef char *charP;
  mc_tpname= new charP[mc_np];

  if(!mc_tpx || !mc_tpname)serror("Error allocating memory\n");
  for(i=0;i<mc_np;i++){
    mc_tpname[i]=new char[20];
    if(!mc_tpname[i])serror("Error allocating memory\n");
  }
  
  if(mc_type==ON_SPACE){
    if(strstr(tmpstr,"random")){
      mc_rnd=1;
      for(i=0;i<mc_np;i++){
        sprintf(mc_tpname[i],"point%d",i);
        mc_tpx[i]=Vector_3(anL*random2()/2,anL*random2()/2,anL*random2()/2);
      }
    }
    else if(strstr(tmpstr,"coords")){
      Set_position("Test coords:");
      for(i=0;i<mc_np;i++){
        Read_param("$: %lf %lf %lf",&mc_tpx[i][0],&mc_tpx[i][1],&mc_tpx[i][2]);
        strncpy(mc_tpname[i],AcFormat,20);
      }
    }
    else serror("Unrecognized specification %s\n",tmpstr);
    
    Set_stop(0);
    mc_ntype=0;
    if(Read_param("Space_particle: %s"),tmpstr){
      if(strstr(tmpstr,"electron"))mc_ntype=0x1;
    }
    Set_stop(1);
      
	    
  }
  else if(mc_type==ON_ION){
    for(i=0;i<mc_np;i++)sprintf(mc_tpname[i],"ion%d",i);
  }
  else if(mc_type==ON_ELC){
    for(i=0;i<mc_np;i++)sprintf(mc_tpname[i],"electron%d",i);
  }

  Read_param("Distribution range: %lf, %lf, %d",&mc_r0, &mc_r1, &mc_nr);
  mc_F0=2.603*pow(npart/(anL*anL*anL),2./3.); // base microfield
  
  Read_param("Print results for points: %s", tmpstr);
  if(mc_plist)delete mc_plist;
  if(strstr(tmpstr,"no"))mc_plist=new /*(cList &)*/cList(0,0);
  else if(strstr(tmpstr,"all"))mc_plist=new /*(cList &)*/cList(0,mc_np);
  else mc_plist=new /*(cList &)*/cList(tmpstr);

  //Read_param("Print average for all points: %s", tmpstr);
  //if(strstr(tmpstr,"y"))mc_prav=1;

  Read_param("Print coord distributions: %s",tmpstr);
  if(strstr(tmpstr,"y"))mc_prcoord=1;

  Read_param("Create files for gnuplot: %s",tmpstr);
  if(strstr(tmpstr,"y"))mc_xmgr=1;

  Set_stop(0);
  if(Read_param("Microfield autocorrelation: %s",tmpstr)){
    if(strstr(tmpstr,"y")){
      printf("Microfield: fourier degree: %d",(int)log2(nsteps)+1);
      ac_n=1<<((int)log2(nsteps)+1);
      printf(", %d points\n",ac_n);
      mc_correl=1;     
    }
  }
  
  mc_force=0x1;
  if(Read_param("Force distribution:>",tmpstr)){
    if(strstr(tmpstr,"coulomb"))mc_force|=COULOMB;
    if(strstr(tmpstr,"effect"))mc_force|=EFFECTIVE;
    if(strstr(tmpstr,"all"))mc_force|=ALL;
  }
  Set_stop(1);

  if(mc_force&EFFECTIVE){
     Read_param("Potential: %s",tmpstr);
     if(strstr(tmpstr,"Kelbg"))mc_potential=PotentialKELBG;
     else if(strstr(tmpstr,"Lennard-Johhnes"))mc_potential=PotentialJONES;
     else if(strstr(tmpstr,"Erf"))mc_potential=PotentialERF;
     else if(strstr(tmpstr,"Cutoff")){
       mc_potential=PotentialCUT;
       Read_param("Cutoff value*: %lf",&E_cut);
     }
     else if(strstr(tmpstr,"ln")){
       mc_potential=PotentialLN;
       Read_param("Cutoff value*: %lf",&E_cut);
     }
     else if(strstr(tmpstr,"table")){
       mc_potential=PotentialTAB;
       Read_param("Potential table file: %s",tmpstr);
       char cfgfile[250];
       strcpy(cfgfile,Get_cur_filename());
       Close_param_file();
       ReadPotential(tmpstr);
       Open_param_file(cfgfile);
     }
     else serror("Unknown potential type specified!\n");


     int pot_corr=0;
     Read_param("Potential correction: %s",tmpstr);
     if(strstr(tmpstr,"y"))
       pot_corr=0x7;
     else{
       if(strstr(tmpstr,"ei") || strstr(tmpstr,"ie"))
         pot_corr|=0x4;
       if(strstr(tmpstr,"ee"))
         pot_corr|=0x1;
       if(strstr(tmpstr,"ii"))
         pot_corr|=0x2;
     }
     if(pot_corr)
       printf("The potential will be corrected in %s%s%s%s%s parts.\n",
         (pot_corr&0x1 ? "ee" : ""),
         (pot_corr&0x2 ? (pot_corr&0x1 ? ", ii" : "ii"): ""),
         (pot_corr&0x4 ? (pot_corr&0x3 ? ", ei" : "ei"): "") );
     else
       printf("The potential will not be corrected.\n");
  
  
     double Clam_ee=1.,Clam_ep=1.;
     
     double Lambda=0.179*sqrt(anT);
     equal_lambda=1;
  
     Lambda_pauli=0.179*sqrt(anT);
     //non_symm=0;
     Lambda_ee=Lambda_pp=Lambda*Clam_ee;
     if(fabs(Clam_ee-1.)>1e-5)equal_lambda=0;
     Lambda_ep=Lambda*Clam_ep;

     if(pot_corr){
       Correction(IONION,mc_potential,anT);
       Correction(IONELC,mc_potential,anT);
       Correction(ELCELC,mc_potential,anT);
     }


  }

  if(mc_force&COULOMB){
    Mf=new Distrvect[mc_np];
    if(!Mf)serror("Microfield_init: memory allocation error.\n");
  }
  if(mc_force&EFFECTIVE){
    Ef=new Distrvect[mc_np];
    //Ef=new Distrvect[mc_np];
    if(!Ef)serror("Microfield_init: memory allocation error.\n");
  }


  for(j=0;j<mc_np;j++){
    for(i=0;i<3;i++){
      if(mc_force&COULOMB)Mf[j][i].init(mc_r0,mc_r1,mc_nr);
      if(mc_force&EFFECTIVE)Ef[j][i].init(mc_r0,mc_r1,mc_nr);
    }
    if(mc_force&COULOMB)Mf[j][3].init(0,mc_r1*sqrt(3.),mc_nr);
    if(mc_force&EFFECTIVE)Ef[j][3].init(0,mc_r1*sqrt(3.),mc_nr);
  }

  if(mc_correl){
    if(mc_force&COULOMB){
      MfCorr=new Correlation[3*mc_np];
      if(!MfCorr)serror("Microfield_init: memory allocation error.\n");
    }
    if(mc_force&EFFECTIVE){
      EfCorr=new Correlation[3*mc_np];
      if(!EfCorr)serror("Microfield_init: memory allocation error.\n");
    }

    for(i=0;i<mc_np;i++){
      for(j=0;j<3;j++){
        if(mc_force&COULOMB){
          MfCorr[3*i+j].init(ac_n,nsteps*dt);
          MfCorr[3*i+j].insert_begin(nsteps);
        }
        if(mc_force&EFFECTIVE){
          EfCorr[3*i+j].init(ac_n,nsteps*dt);
          EfCorr[3*i+j].insert_begin(nsteps);
        }
      }
    }
  }
}





void Microfield_step(){
  int i,j;
  Vector_3 fe(0,0,0), fc(0,0,0),r;
  double R,f;
  double dEpotent, dQuant;

  if(mc_type==ON_ELC)mc_tpx=anx+nion;
  if(mc_type==ON_ION)mc_tpx=anx;

  if(!xx_used){
    xx_used=1;
    for(i=0;i<npart;i++)anx[i]=r_cell(anx[i]);
  }

  
  int type=0;
  for(j=0;j<mc_np;j++){

    if(mc_type!=ON_SPACE){
      if(j>=nion)type|=0x1;// setting electron-? interaction type
      type&=0x1;         // setting ?-ion interaction type
    }
    else type=mc_ntype;

  
    fe=Vector_3(0,0,0);
    fc=Vector_3(0,0,0);
    for(i=0;i<npart;i++){

      if(i>=nion)type|=0x2;// setting ?-electron interaction type 
      
      if(mc_type==ON_ELC && i==j+nion)continue;
      if(mc_type==ON_ION && i==j)continue;
       
      r=anx[i]-mc_tpx[j];
      
      int k;
      for(k=0;k<3;k++){
	if(r[k]>anL/2)r[k]-=anL;
	if(r[k]<-anL/2)r[k]+=anL;
      }

    
     
      //r.print();
      R=r.normalize();
      if(R>1e-10){
	if(mc_force&COULOMB){
	  double a=R*R*mc_F0;
	  if(i<nelec)fc+=r/a;
	  else fc-=r/a;
	}
	if(mc_force&EFFECTIVE){
	  f=mc_potential(type,R,dEpotent,dQuant);
	  fe+=r*f/mc_F0;
	}
      }
      else printf("(%d %d)!\n",j,i);
    }

    if(mc_force&COULOMB){
      for(i=0;i<3;i++)Mf[j][i].point(fc[i],1.);
      Mf[j][3].point(fc.norm(),1.);
    }
    if(mc_force&EFFECTIVE){
      for(i=0;i<3;i++)Ef[j][i].point(fe[i],1.);
      Ef[j][3].point(fe.norm(),1.);
    }

    if(mc_correl){
      if(mc_force&COULOMB){
	      for(i=0;i<3;i++)MfCorr[3*j+i].insert_next(fc[i]);
      }
      if(mc_force&EFFECTIVE){
	      for(i=0;i<3;i++)EfCorr[3*j+i].insert_next(fe[i]);
      }
    }

  }
}   

# if 0

double Holtsmark(double x){
  if(x<=1e-6)return 0.;

  double dt=0.01*2*M_PI/x;
  double t=dt,sum=0.,ds;
  do{
   ds=exp(-pow(t,3./2.))*sin(x*t);
   sum+=ds;
  }while(fabs(ds)>1e-5);

  return (2/M_PI)*x*sum;
}

# endif

void output_mfdistr(char *filename, char *dname, Distrvect *Mf){
  char str[250];
  int i;
  FILE *f1=Err_fopen(filename,"wt"), *f2=NULL;

  if(mc_xmgr){
    strcpy(str,out_dir);
    strcat(str,dname);
    strcat(str,"-mf.gnu");
    f2=Err_fopen(str,"wt");
  }
  
    // creating header
  fprintf(f1,"#// Microfield distribution data\n");
  fprintf(f1,"#Test points (coords in cell length L):\n");
  for(i=0;i<mc_np;i++){
    fprintf(f1,"#%s: %f, %f, %f\n",mc_tpname[i],mc_tpx[i][0]/anL,mc_tpx[i][1]/anL,mc_tpx[i][2]/anL);
  }
  fprintf(f1,"\n\n#%13s","   ");
  
  int j;
  int lfield=13;
  if(mc_prcoord)lfield=13*4;
  char frm[20];
  sprintf(frm,"%%%ds",lfield);
  while((j=mc_plist->step())>=0){
    fprintf(f1,frm,mc_tpname[j]);
  }
  fprintf(f1,frm,"average");
  fprintf(f1,frm,"average_c");
  //fprintf(f1,frm,"Holtsmark");
  fprintf(f1,"\n#%13s","Fx (Fr)");

  
  char xyz[]="Wx\0Wy\0Wz\0r";
  for(i=0;i< mc_plist->n + 1;i++){
    if(mc_prcoord)for(i=0;i<3;i++)fprintf(f1,"%13s",xyz+3*i);
    fprintf(f1,"%13s","Wr*F^2");
  }
  fprintf(f1,"%13s","Wr*F^2");
  fprintf(f1,"%13s","Wr*F^2");
  //fprintf(f1,"%13s","Wr*F^2");
  fprintf(f1,"\n");
  
  Statistics st[5];
 
 
  double dx, lim=mc_r1*sqrt(3.);
  double x, x0=(mc_prcoord ? mc_r0 : 0);

  /*   No normalization !!!
  for(i=0;i<mc_np;i++){
    for(j=0;j<4;j++)Mf[i][j].normalize(1.);
    //Mf[i][2].normalize_r2(1/(4*M_PI));
  } */
  

  dx=(lim-x0)/mc_nr/5;
  
  for(x=x0;x<=lim;x+=dx){

   for(j=0;j<5;j++)st[j].clear();
   for(i=0;i<mc_np;i++){
     for(j=0;j<3;j++){
      st[j].next((Mf[i][j])(x),1);
      st[4].next(Mf[i][j].nrm_r2(x)+Mf[i][j].nrm_r2(-x),1);
     }
     st[3].next((Mf[i][3])(x),1);
   }
  
   fprintf(f1,"%12f ",x);
   if(mc_xmgr)fprintf(f2,"%12f ",x);

   mc_plist->rewind();
   while((j=mc_plist->step())>=0){
     if(mc_prcoord){
       for(i=0;i<3;i++){
	 fprintf(f1,"%12f ",(Mf[j][i])(x));
	 if(mc_xmgr)fprintf(f2,"%12f ",(Mf[j][i])(x));
       }
       //fprintf(f1,"%12f ",(Mf[2][i])(x)*x*x);
       //if(mc_xmgr)fprintf(f2,"%12f ",(Mf[2][i])(x)*x*x);
     }
     fprintf(f1,"%12f ",(Mf[j][3])(x));
     if(mc_xmgr)fprintf(f2,"%12f ",(Mf[j][3])(x));
   }
   if(mc_prcoord){
     for(i=0;i<3;i++){
       fprintf(f1,"%12f ",st[i].av());
       if(mc_xmgr)fprintf(f2,"%12f ",(double)st[i].av());
     }
     //fprintf(f1,"%12f ",st[2].av()*x*x);
     //if(mc_xmgr)fprintf(f2,"%12f ",st[2].av()*x*x);
   }
   fprintf(f1,"%12f",st[3].av());
   if(mc_xmgr)fprintf(f2,"%12f",(double)st[3].av());
   fprintf(f1,"%12f\n",st[4].av());
   if(mc_xmgr)fprintf(f2,"%12f\n",(double)st[4].av());   
   //fprintf(f1,"%12f  \n",Holtsmark(x));
   //if(mc_xmgr)fprintf(f2,"%12f \n",Holtsmark(x));
  }

  if(mc_xmgr)fclose(f2); 
  fclose(f1);
}

void output_mfcorr(char *filename,Correlation* MfCorr){
  int i;
  FILE *f1=Err_fopen(filename,"wt");
  printf("Microfiel corelations: performing  FFT...\n");
  double *fft_sq=Calculate_forth(MfCorr[0].n,MfCorr,3*mc_np);
  printf("Inverse FFT...\n");
  Calculate_back(fft_sq,MfCorr[0].n);
  for(i=0;i<MfCorr[0].n;i++)fft_sq[i]=fft_sq[2*i];
  TableFunction acf(MfCorr[0].n,fft_sq);
  double tmax=MfCorr[0].maxtime*Wp;
  acf.xscale(0.,tmax);
  double dt=0.05,t;
  for(t=0.;t<tmax && t< 150.;t+=dt){
    fprintf(f1,"%f  %f\n",t,acf(t)/acf(0.));
  }
  fclose(f1);
}

void Microfield_process(char *dname){
  char str[256]="microf.dat";
 
  if(mc_force&COULOMB){
    strcpy(str,out_dir);
    strcat(str,dname);
    if     (mc_type==ON_ION)strcat(str,"-mfi.dat");
    else if(mc_type==ON_ELC)strcat(str,"-mfe.dat");
    else strcat(str,"-mfn.dat");
    output_mfdistr(str,dname, Mf);
    
    if(mc_correl){
      strcpy(str,out_dir);
      strcat(str,dname);
      if     (mc_type==ON_ION)strcat(str,"-EEi.dat");
      else if(mc_type==ON_ELC)strcat(str,"-EEe.dat");
      else strcat(str,"-EEn.dat");
      output_mfcorr(str,MfCorr);   
    }
  }  
  if(mc_force&EFFECTIVE){
    strcpy(str,out_dir);
    strcat(str,dname);
    if     (mc_type==ON_ION)strcat(str,"-efi.dat");
    else if(mc_type==ON_ELC)strcat(str,"-efe.dat");
    else strcat(str,"-efn.dat");
    output_mfdistr(str,dname, Ef);
    
    if(mc_correl){
      strcpy(str,out_dir);
      strcat(str,dname);
      if     (mc_type==ON_ION)strcat(str,"-ffi.dat");
      else if(mc_type==ON_ELC)strcat(str,"-ffe.dat");
      else strcat(str,"-ffn.dat");
      output_mfcorr(str,EfCorr);   
    }
  } 
}  
 



double cur_tm=-1;
long cur_n;
char cur_ffile[200]="";
FILE *cur_ff1=NULL;
double *cur_arr;
long cur_k;
Vector_3 *cur_vec;
Statistics *cur_sum;


# define ccc(i,j) cur_arr[(long)(2*cur_n-(i)-1)*(long)(i)/2+(j)]





int cur_f=0;
int cur_comp=0;
int cur_part=0;

Vector_3 *vold;
Vector_3 *jmean;

cList* cur_vv_list=NULL;
int cur_avr;

void Cur_init(){
  if(!(tType&FLOW) && !(tType&VEL))serror("Can't calculate current-current "
   "correlations,\n the trajectory file does not contain appropriate data.\n");

  Set_stop(0);
  Read_param("File for current data: %s",cur_ffile);
  Read_param("Max. time: %lf",&cur_tm);
  if(cur_tm<0 || cur_tm>=dt*nsteps ){
      printf("Warning: maximal time %f is greater than trajectory length %f\n"
             "or not specified, setting to full length\n",cur_tm,dt*nsteps);
      cur_tm=nsteps*dt;
  }
  cur_n=(long)(cur_tm/dt);


  Set_stop(1);
  char tmpstr[200];

  cur_part=1;
  Read_param("Components: %s",tmpstr);
  if(strstr(tmpstr,"all"))cur_comp=0;
  else if(strstr(tmpstr,"e-e"))cur_comp=1;
  else if(strstr(tmpstr,"e-p"))cur_comp=1;
  else{
   if(cur_vv_list)delete cur_vv_list;
   cur_vv_list=new /*(cList &)*/cList(tmpstr);
   if(!cur_vv_list->valid())serror("Undefined current component!\n");
   cur_vv_list->adjust(0,npart-1);
   cur_comp=-1;
   cur_vv_list->rewind();
   cur_part=cur_vv_list->n;
   if(cur_part<=0)serror("Void components list !\n");
   if(cur_part>=npart)serror("Error adjusting components list!\n");
  }

  if(!time_aver){ // no time averaging
    cur_arr= new double[nsteps*cur_part];
    vold= new Vector_3[cur_part];

    if(!cur_arr || !vold )serror("Cur_init: Memory allocation error!\n");
    cur_k=0;

    return;
  }

  Read_param("Processing type: %s",tmpstr);

  if( strstr(tmpstr,"four") ){
    cur_f=1;

    //Read_param("Fourier degree: %ld",&cur_n);
    printf("Fourier degree: %d",(int)log2(cur_n)+1);
    cur_n=1<<((int)log2(cur_n)+1);
    printf(", %ld points\n",cur_n);

    cur_arr= new double[6*cur_part*cur_n];
    jmean = new Vector_3[cur_part];
    vold=   new Vector_3[cur_part];

    if(!cur_arr || !jmean || !vold)serror("Cur_init: Memory allocation error!\n");

    cur_k=0;
    long i;
    for(i=0;i<6*cur_n*cur_part;i++)cur_arr[i]=0;

  }
  else if(strstr(tmpstr, "correct")){
    cur_f=2;
    cur_arr= new double[3*nsteps];
    jmean = new Vector_3;
    vold=   new Vector_3;
    if(!cur_arr)serror("Cur_init: Memory allocation error!\n");
    cur_k=0;
  }
  else{
    cur_f=0;

    //Set_stop(0);

    //Set_stop(1);
    //Read_param("Time intervals: %d",&cur_nt);

    cur_arr= new double[cur_n*(cur_n+1)/2];
    cur_sum= new Statistics[cur_n];
    cur_vec= new Vector_3[cur_n];
    jmean = new Vector_3;
    vold=   new Vector_3;
    if(!cur_arr || !cur_sum || !cur_vec)serror("Cur_init: Memory allocation error!\n");

    cur_k=cur_n;
  }
  Set_stop(0);
  cur_avr=1;
  if(Read_param("Relative to average: %s",tmpstr)){
   if(strstr(tmpstr,"no"))cur_avr=0;
  }


  if(cur_ffile[0])cur_ff1=Err_fopen(cur_ffile,"wt");

  for(int i=0;i<cur_part;i++){
    for(int j=0; j<3;j++)jmean[i][j]=vold[i][j]=0.;
  }
}


double stpx=0;
int stpi=0;


void add_stp_notime(Vector_3 &vs){
  if(stpi==0){
    vold[0]=vs;
  }
  cur_arr[stpi]=(vold[0]*vs)/(vold[0]*vold[0]);
  stpi++;
}

void add_stpvv_notime(){
  int i;
  cur_vv_list->rewind();
  if(stpi==0){
     for(i=0;i<cur_part;i++){
      int ipart=cur_vv_list->step();
      vold[i]=anv[ipart];
     }
     cur_vv_list->rewind();
  }
  cur_arr[stpi]=0;
  for(i=0;i<cur_part;i++){
   int ipart=cur_vv_list->step();
   cur_arr[stpi]+=(vold[i]*anv[ipart])/(vold[i]*vold[i]);
  }
  cur_arr[stpi]/=cur_part;
  stpi++;
}



void add_stp(Vector_3 &vs){

  if(!time_aver){
    add_stp_notime(vs);
    return;
  }

  static Vector_3 vr;

  if(cur_f==1){
   double d;
   int dumm;


   while(stpi<cur_n && (double)cur_k>stpx){

    d=DxScale(1.,stpx,dumm);
    vr=vold[0]*d+vs*(1-d);
    jmean[0]+=vr;

    cur_arr[2*stpi]=        vr[0];
    cur_arr[2*cur_n+2*stpi]=vr[1];
    cur_arr[4*cur_n+2*stpi]=vr[2];

    stpi++;
    stpx=((double)(nsteps-1))*stpi/(cur_n/2-1);

   }
   vold[0]=vs;
   cur_k++;
   return;
  }
  if(cur_f==2){
   jmean[0]+=vs;
   cur_arr[cur_k]=vs[0];
   cur_arr[nsteps+cur_k]=vs[1];
   cur_arr[2*nsteps+cur_k]=vs[2];

   cur_k++;
   return;
  }

  long i,j;
  for(i=cur_k;i<cur_n;i++){
    for(j=i;j<cur_n;j++){
     if(i==0 && cur_k==0)cur_sum[j].next(ccc(0,j));
     else ccc(i-1,j-1)=ccc(i,j);
    }
    if(i!=0){
      cur_vec[i-1]=cur_vec[i];
      ccc(i-1,cur_n-1)=vs*cur_vec[i];
    }
  }
  cur_vec[cur_n-1]=vs;
  ccc(cur_n-1,cur_n-1)=vs*vs;

  if(cur_k>0)cur_k--;

}

void add_stpvv(){
  if(!time_aver){
    add_stpvv_notime();
    return;
  }
  if(cur_f!=1)serror("Can make only fourrier <vv> calculations!\n");
  long i;
  double *tmparr=cur_arr;

  Vector_3 vr;
  double d;
  int dumm;

  while(stpi<cur_n && (double)cur_k>stpx){

    d=DxScale(1.,stpx,dumm);
    cur_vv_list->rewind();
    for(i=0;i<cur_part;i++){
      int ipart=cur_vv_list->step();

      cur_arr=tmparr+6*cur_n*i;

      vr=vold[i]*d+anv[ipart]*(1-d);

      jmean[i]+=vr;

      cur_arr[2*stpi]=        vr[0];
      cur_arr[2*cur_n+2*stpi]=vr[1];
      cur_arr[4*cur_n+2*stpi]=vr[2];

    }

    stpi++;
    stpx=((double)(nsteps-1))*stpi/(cur_n/2-1);


  }
  cur_vv_list->rewind();
  for(i=0;i<cur_part;i++){
   int ipart=cur_vv_list->step();
   vold[i]=anv[ipart];
  }
  cur_k++;

  cur_arr=tmparr;
  return;
}




void Flow(int n, Vector_3 *v, Vector_3* flowe, Vector_3* flowp){
  int nn2=n/2;
  int i;
  for(i=0;i<nn2;i++)*flowe+=v[i];

  for(i=nn2;i<n;i++)*flowp+=v[i];

  (*flowe)/=nn2;
  (*flowp)/=nn2;
}



void Cur_step(){
//int  nn2;
//nn2=npart/2;

  Vector_3 vs(0,0,0);

  if((tType&VEL) && cur_comp>=0)Flow(npart,anv,flow,flow+1);

  if(cur_comp>0)vs=(flow[0]+flow[1])/2.;
  else if(cur_comp<0){
    add_stpvv();
    return;
  }
  else vs=flow[0]-flow[1];

  add_stp(vs);
  if(cur_ff1)fprintf(cur_ff1,"%13e %13e %13e %13e %13e\n",cur_step*dt*Wp,vs[0],vs[1], vs[2], vs*vs);

}

void calc_corr(double *arr,long lag,long n,double *corr,Vector_3& jmean){
  long i,k;

  printf("%ld, %ld\n", n,lag);
//  for(k=0;k<3;k++){
//   jmean[k]=0;
//   for(i=0;i<n;i++)/* printf("%f-",arr[k*n+i]),*/jmean[k]+=arr[k*nsteps+i];
//   jmean[k]/=n;
//  }
  jmean/=n;
  int pr=lag/100;
  for(i=0;i<lag;i++){
   if(!(i%pr))printf("%.0f %%\n",(double)i*100./lag);
   corr[i]=0;
   double cmean=0, csq=0;
   for(k=0;k<n-i;k++){
    corr[i]+=(arr[k]-jmean[0])*(arr[k+i]-jmean[0]);
    corr[i]+=(arr[n+k]-jmean[1])*(arr[n+k+i]-jmean[1]);
    corr[i]+=(arr[2*n+k]-jmean[2])*(arr[2*n+k+i]-jmean[2]);
    cmean+=corr[i];
    csq+=corr[i]*corr[i];
   }
   corr[lag+i]=sqrt((csq-cmean*cmean/(n-i))/(n-i));
  }
  double var=corr[0];
  for(i=0;i<lag;i++){
   corr[i]/=var*(n-i)/n;
   corr[lag+i]/=var*(n-i)/n;
  }
}

extern double start_time;

void Cur_process(char *dname){
  if(cur_ff1)fclose(cur_ff1);
  printf("Calculating j-j correlation...\n");
  long i,l;

  strcpy(cur_ffile,out_dir);
  strcat(cur_ffile,dname);
  strcat(cur_ffile,"-jj.dat");


  FILE *f1=Err_fopen(cur_ffile,"wt");

  if(!time_aver){
   printf("Writing data...\n");
   for(i=0;i<nsteps;i++){
     fprintf(f1,"%13e   %13e\n",Wp*dt*i,cur_arr[i]);
   }
   fclose(f1);
   return;
  }

  if(!cur_avr){
   for(l=0;l<cur_part;l++){
    for(i=0;i<3;i++)jmean[l][i]=0;
   }
  }

  if(cur_f==1){
    double *tmp=cur_arr;
    if(cur_avr){
     for(l=0;l<cur_part;l++){
      cur_arr=tmp+6*cur_n*l;
      int j;
      for(i=0;i<3;i++)jmean[l][i]=0;
      for(j=0;j<cur_n/2;j++){
       for(i=0;i<3;i++)jmean[l][i]+=cur_arr[2*cur_n*i+2*j];
      }
     }
    }



    for(l=0;l<cur_part;l++){
      printf("Paricle %ld:\n",l);
      cur_arr=tmp+6*cur_n*l;
      jmean[l]/=(cur_n);     // !!!!!!

      for(i=0;i<3;i++){
       int j;
       for(j=0;j<cur_n/2;j++){
         cur_arr[2*cur_n*i+2*j]-=jmean[l][i];
         cur_arr[2*cur_n*i+2*(cur_n-j-1)]= cur_arr[2*cur_n*i+2*j];
         cur_arr[2*cur_n*i+2*j+1]=cur_arr[2*cur_n*i+2*(cur_n-j-1)+1]=0.;
       }
       printf("Fourier: component #%ld\n",i+1);
       four(cur_arr+2*i*cur_n,cur_n,1);
      }

  				for(i=0;i<cur_n;i++){


       cur_arr[2*i]=cur_arr[2*i]*cur_arr[2*i]+cur_arr[2*i+1]*cur_arr[2*i+1];

       cur_arr[2*i]+=cur_arr[2*(i+cur_n)]*cur_arr[2*(i+cur_n)]+cur_arr[2*(i+cur_n)+1]*cur_arr[2*(i+cur_n)+1];
       cur_arr[2*i]+=cur_arr[2*(i+2*cur_n)]*cur_arr[2*(i+2*cur_n)]+cur_arr[2*(i+2*cur_n)+1]*cur_arr[2*(i+2*cur_n)+1];

       cur_arr[2*i+1]=0;
       //cur_arr[2*(cur_n+i)]=cur_arr[2*i]; //copying fourier transform square



      }
    }
    cur_arr=tmp;
    for(l=1;l<cur_part;l++){
      cur_arr=tmp+6*cur_n*l;
      for(i=0;i<cur_n;i++){
       tmp[2*i]+=cur_arr[2*i];
      }
    }

    for(i=0;i<cur_n;i++){
      tmp[2*i]/=cur_part;
      tmp[2*i+1]=0.;
      tmp[2*(cur_n+i)]=tmp[2*i];//copying fourier transform square
    }
    cur_arr=tmp;


    printf("Back fourier\n");

    four(cur_arr,cur_n,-1);
    printf("Writing data...\n");
    for(i=0;i<cur_n/2;i++){
      fprintf(f1,"%13e   %13e  %13e  %13e\n",Wp*dt*i*(nsteps-1)/(cur_n/2-1),cur_arr[2*i]/cur_arr[0],2*M_PI*(double)i/(Wp*dt*(nsteps-1)),cur_arr[2*(i+cur_n)]);
    }
    fclose(f1);
    return;
  }


  if(cur_f==2){
    cur_n=(cur_step < cur_n ? cur_step : cur_n);
    double *corr=new double[2*cur_n];
    if(!corr)serror("Error allocating memory\n");
//    Vector_3 jmean(0,0,0);
    printf("Processing...\n");
    calc_corr(cur_arr,cur_n,cur_step,corr,jmean[0]);
    printf("Writing data...\n");
    for(i=0;i<cur_n;i++){
      fprintf(f1,"%13e   %13e   %13e\n",Wp*dt*i,corr[i],corr[cur_n+i]);
    }
   //fprintf(f1, "Jmean: -  -  - %f, %f, %f\n",jmean[0],jmean[1],jmean[2]);
    fclose(f1);
    delete corr;
    return;
  }

  Vector_3 vv(0,0,0);
  for(i=0;i<cur_n;i++){
    add_stp(vv);
  }

  double c0=cur_sum[0].s, t=start_time;
  cur_sum[0]/=c0;
  //printf("%ld -------\n",nsteps);
  double er0=cur_sum[0].rel(cur_step)/sqrt((double)cur_step);
  double eri;
  //printf("%f\n",c0);
  for(i=0;i<cur_n;i++){
    if(i>0){
     cur_sum[i]/=c0;
     eri=cur_sum[i].rel(cur_step-i)/sqrt((double)cur_step-i)+er0;
    }
    else eri=er0;
    double corr=cur_sum[i].s;
    fprintf(f1,"%13e   %13e  %13e\n",t*Wp,corr,eri*fabs(corr));
    t+=dt;
  }
  fclose(f1);
}


FILE *new_trj=NULL;
long ref_type=0;
extern double max_time;

void Refine_init(char *idname){
  char tmpstr[200],str[200];
  Set_stop(1);
  Read_param("New trajectory file: %s",tmpstr);
  sprintf(str,tmpstr,idname); // setting new file name
  strcpy(tmpstr,str);
  strcpy(str,out_dir);
  strcat(str,tmpstr);
  new_trj=Err_fopen(str,"wb");

  Set_stop(0);
  Read_param("Write data on:>",tmpstr);


  if(strstr(tmpstr,"vel")){
    if(!(tType&VEL))printf("The trajectory file does not contain velocities.\n");
    else ref_type|=VEL;
  }
  if(strstr(tmpstr,"coord")){
    if(!(tType&COORD))printf("The trajectory file does not contain coordinates.\n");
    else ref_type|=COORD;
  }

  if(strstr(tmpstr,"flow"))ref_type|=FLOW;

  if(max_time>0. && !(tType &COORD&VEL))
    printf("Refine: warning: max_time is set!!!\n");


  fseek(tFile,0,SEEK_SET);
  char buffer[200];
  fread (buffer,sizeof(char),200,tFile);
  fwrite(buffer,sizeof(char),200,new_trj);

  fread(buffer,sizeof(char),200,tFile); // reading

  int nn;
  double l,g,t,dt0,Wp,mass=1.;
  long nsteps0, type;
  char dname[100];
  sscanf(buffer,"%d %lf %lf %lf %s %lf %ld %lf %ld %lf", &nn, &l,&g,&t,dname,&dt0,&nsteps0,&Wp,&type,&mass);
  type=ref_type;
  long newnstp=nsteps0;
  if(fabs(dt0/dt-1.)>1e-5){
    newnstp=(long)(nsteps0/dt*dt0);
    printf("Refine: adjusting dt X nsteps...\n"
	   "old:    %f X %ld =%f\n"
	   "new:    %f X %ld =%f\n",
	   dt0,nsteps0,dt0*nsteps0,dt,newnstp,dt*newnstp);
  }
  sprintf(buffer,"%d %f %f %f %s %f %ld %f %ld %f", nn, l,g,t,dname,dt,newnstp,Wp,type,mass);
  fseek(tFile,400,SEEK_SET);
  fwrite(buffer,sizeof(char),200,new_trj);
}


void Refine_step(){
  if(ref_type&FLOW){
   Flow(npart,anv,flow,flow+1);
   fwrite(flow,sizeof(Vector_3),2,new_trj);
  }
  int nw;
  //long pos;
  if(ref_type&COORD){

    nw=fwrite(anx,sizeof(Vector_3),npart,new_trj);
    if(nw!=npart)printf("Write failed!\n");

  }
  if(ref_type&VEL){
    nw=fwrite(anv,sizeof(Vector_3),npart,new_trj);
    if(nw!=npart)printf("Write failed!\n");

  }

}

void Refine_process(){
  fclose(new_trj);
}


# define MAXLEV 5
# define MAXSPLITS 200


unsigned arrp=0;

int _splitsq(unsigned n, int lev, unsigned *arr,unsigned limit){
  double lsq=n;
  lsq=sqrt(lsq);  
  unsigned rr,rest;
  rr=(unsigned)lsq;

  rest=n-rr*rr;
  if(lev==0 && rest !=0)return 0;

  int nsp=0,i;
  do{
    arr[arrp*MAXLEV+lev]=rr;
    if(rest==0){
      for(i=lev-1;i>=0;i--)arr[arrp*MAXLEV+i]=0;
      if(arrp>=limit-1)printf("splitsq: Array limit reached!\n");
      else arrp++;
      for(i=MAXLEV-1;i>=lev;i--)arr[arrp*MAXLEV+i]=arr[(arrp-1)*MAXLEV+i];
      nsp++;
      if(lev==0 || rr==0)return nsp;
    }
    else nsp+=_splitsq(rest,lev-1,arr,limit);

    if(rr==0)break;
    rr--;
    rest=n-rr*rr;

  }while(1/*rr!=0 && rest<(rr+1)*(rr+1)*/);

  return nsp;
} 




int splitsq(unsigned n, int lev,unsigned *arr,unsigned limit=MAXSPLITS){
  arrp=0;
  if(lev>MAXLEV){
    printf("spltsq: Level is reduced to MAXLEV=%d\n",MAXLEV);
    lev=MAXLEV;
  }
  return _splitsq(n,lev-1,arr,limit);
}    



  
int statstruc=0;
Statistics *R0ee, *R0ep;
double ks,ke;
Vector_3 kd;
int kn,kvps=3, kdir=0;
int *kvec_list=NULL;

void Sstruc_init(){
  char tmpstr[100];
  if(!(tType&COORD))serror("Can't calculate static structure factor,\n"
			   " trajectory file does not contain Coords data.\n");
 
  Set_stop(0);
  if(Read_param("K direction :>",tmpstr)){
    if(strstr(tmpstr,"allowed")){
      kdir=2;
      Set_stop(1);     
      Read_paramn(1,"Vectors per step: %d",&kvps);
      
    }  
    else if(sscanf(tmpstr,"%lf,%lf, %lf",&kd[0],&kd[1],&kd[2])==3){
      kdir=1;
      kd.normalize();
      kvps=1;
    }
    else if(!strstr(tmpstr,"ort"))serror("Undefined K direction.\n");
  }
  else {
    kdir=0;
    kvps=3;
    printf("No K-direction specified, using ort-averaging.\n");
  }
  Set_stop(1);
  Read_paramn(3,"K interval, steps: %lf,%lf,%d",&ks,&ke,&kn);

  if(kdir==2){ // preparing vector list

    int i;
    kvec_list= new int[3*kn*kvps];
    if(!kvec_list)serror("Error allocating k-vector list.\n");
    for(i=0;i<3*kn*kvps;i++){
      kvec_list[i]=0;
    }
      

    double k, kstp;
    if(kn>=1)kstp=(ke-ks)/(kn-1);
    else kstp=0;

    for(i=0,k=ks;i<kn;i++,k+=kstp){
      unsigned v1=(unsigned)(k*k)+1;
      unsigned v2=(unsigned)((k+kstp)*(k+kstp));
      unsigned nv=v2-v1+1;

      if(nv==0)continue;
      unsigned *arr=new unsigned[MAXLEV*MAXSPLITS*nv];
      if(!arr)serror("Structure factor vector preparing: MAE\n");
      
      unsigned j;
      unsigned l=0;
      for(j=v1;j<=v2;j++){
	l+=splitsq(j,3,arr+MAXLEV*l);
      }
      int copied=0;
      for(j=0;j<l;j++){ // selecting vectors
	int ac=1;
	if(copied<kvps){
	  if((0.5*(1+random2()))*(l-j)>(kvps-copied))ac=0;
	}
	if(ac){
	  memcpy(kvec_list+i*3*kvps,arr+MAXLEV*j,3*sizeof(int));
	  copied++;
	}
	
      }
      delete arr;
    }
  }
    

  R0ee=new Statistics[kn];
  R0ep=new Statistics[kn];
  if(!R0ee || !R0ep)serror("Can not allocate array for structure factor.\n");
}


int set_k_vector(int i, int j){
  int k;
  int ac=1;
  switch(kdir){
  case 0:
    kd[0]=kd[1]=kd[2]=0;
    kd[j]=1.;
    break;
  case 2:
    ac=0;
    for(k=0;k<3;k++){
      kd[k]=(double)kvec_list[3*kvps*i+3*j+k];
      if(kd[k]!=0.)ac=1;
    }
    break;
  }
  return ac;
}


void step_allowed(){
  double k, kstp;
  int i,j,jj,m,l;
  int nn2=npart/2;
  Vector_3 r;

  if(kn>=1)kstp=(ke-ks)/(kn-1);
  else kstp=0;

  struct crl {
    double dres,drec, drps,drpc;
  } c;

  for(i=0,k=ks;i<kn;i++,k+=kstp){
    //double avee=0,avep=0;

    for(j=0;j<kvps;j++){

      l=set_k_vector(i,j);

      if(l==0){
	//break;
	//avee=avep=0;

      }
      else{
	kd*=2*M_PI/anL;

	int e;
	m=0;
	while(kd[m]==0. && m<3)m++;

	double dee=0, dep=0;
	int aver=0;
	for(e=0;e<4;e++){ // averaging through 8 directions

	  c.drec=c.dres=c.drpc=c.drps=0.;

	  for(jj=0;jj<npart;jj++){
	    if(!xx_used)r=r_cell(anx[jj]);
	    else r=anx[jj];

	    if(jj<nn2){
	      c.drec+=cos(r*kd);
	      c.dres+=sin(r*kd);
	    }
	    else{
	      c.drpc+=cos(r*kd);
	      c.drps+=sin(r*kd);
	    }
	  }

	  dee+=0.5*(c.dres*c.dres+c.drps*c.drps+c.drec*c.drec+c.drpc*c.drpc)/nn2;
	  dep+=(c.drec*c.drpc+c.dres*c.drps)/nn2;

	  aver++;
	  if(e!=0)kd[m]=-kd[m];

	  m=(m+1)%3;
	  if(kd[m]==0.){
	    e+=2;
	    m=(m+1)%3;
	    if(kd[m]==0.)break;
	  }
	  kd[m]=-kd[m];
	}

	dee/=aver;
	dep/=aver;


	R0ee[i].next(dee);
	R0ep[i].next(dep);

      }



    }

  }
}



void Sstruc_step(){
  if(kdir==2){
    step_allowed();
    return;
  }
  int i,j,m;
  int nn2=npart/2;
  double k;
  Vector_3 dres,drec, drps,drpc;
  Vector_3 r;

  for(i=0;i<kn;i++){
    k=(ks+((ke-ks)*i)/(kn-1))*2*M_PI/anL;
    dres=drec=Vector_3(0.,0.,0.);
    for(j=0;j<nn2;j++){
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];
      if(kdir){
	drec[0]+=cos((r*kd)*k);
	dres[0]+=sin((r*kd)*k);
      }
      else{
	for(m=0;m<3;m++){
	  drec[m]+=cos(r[m]*k);
	  dres[m]+=sin(r[m]*k);
	}
      }
    }
    //dre/=srt(nn2);

    drps=drpc=Vector_3(0.,0.,0.);
    for(j=nn2;j<npart;j++){
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];
      if(kdir){
	drpc[0]+=cos((r*kd)*k);
	drps[0]+=sin((r*kd)*k);
      }
      else{
	for(m=0;m<3;m++){
	  drpc[m]+=cos(r[m]*k);
	  drps[m]+=sin(r[m]*k);
	}
      }
    }
    //drp/=nn2;

    double avee=0,avep=0;
    int mm=3;
    if(kdir)mm=1;

    for(m=0;m<mm;m++){
      avee+=0.5*(dres[m]*dres[m]+drps[m]*drps[m]+drec[m]*drec[m]+drpc[m]*drpc[m])/(mm*nn2);
      avep+=(drec[m]*drpc[m]+dres[m]*drps[m])/(mm*nn2);
    }


    R0ee[i].next(avee);
    R0ep[i].next(avep);
  }
}

void Sstruc_process(char *dname){
  char str[256]="statstruc.dat";
  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-ss.dat");
  int i;

  FILE *f1=Err_fopen(str,"wt");
  printf("Writing structure factor...\n");
  double k;
  for(i=0;i<kn;i++){
    k= ks+((ke-ks)*i)/(kn-1);
    fprintf(f1,"%f %f %f %f %f %f\n",
	    k, R0ee[i].av(),R0ep[i].av(),
	    4*M_PI*k*k*(R0ee[i].av()-1.),
	    4*M_PI*k*k*R0ep[i].av(), R0ee[i].av()-R0ep[i].av()-1. );
  }
  fclose(f1);
}




long ds_nf;

typedef double *pfarr;

pfarr *ds_e,*ds_p;
int ds_count=0;
int ds_nvect=0;
Vector_3 *ds_kvect;

double ds_stpx=0;
int ds_stpi=0;

class comple{
public:
  double x;
  double y;
  double &operator[](int i){
    if(i%2==0)return x;
    else return y;
  }
};

comple *eold, *pold, *dre, *drp;
double ds_f0,ds_f1,ds_fc;
int ds_dir=0;

void Dstruc_init(){
  int l,i;

  if(!(tType&COORD))serror("Can't calculate dynamic structure factor,\n"
			   " trajectory file does not contain Coords data.\n");
  unsigned int arr[MAXLEV*MAXSPLITS];

  Set_stop(0);
  
  if(time_aver){
    if(Read_param("Direct Fourier: %lf %lf %lf",&ds_f0,&ds_f1,&ds_fc)==3){
      ds_dir=1;
      printf("Direct transform !\n");
    }

    printf("DSF fourier degree: %d",(int)log2(nsteps)+1);
    ds_nf=1<<((int)log2(nsteps));
    printf(", %ld points\n",ds_nf);
  }
  else{
    ds_nf=nsteps;
  }


  double kval;
  if(Read_param("Free K: %lf %d", &kval, &ds_nvect)==2){
    printf("Warning: using not allowed k-directions! \n");
    ds_kvect= new Vector_3[ds_nvect];
    if(!ds_kvect)serror("DSF: MAE\n");
    for(i=0;i<ds_nvect;i++){
      ds_kvect[i]=Vector_3(random2(),random2(),random2());
      ds_kvect[i].normalize();
      ds_kvect[i]*=kval*2*M_PI/anL;
    }
  }
  else{

    Set_stop(1);
    int ksq;
    Read_param("K square: %d",&ksq);
    //Read_param("Frequency interval, steps: %f, %f, %d",&ds_f0,&ds_f1,&ds_nf);

    int nk;
    nk=splitsq(ksq,3,arr);
    printf("%d components\n",nk);



    ds_kvect= new Vector_3[8*nk];
    if(!ds_kvect)serror("DSF: MAE\n");

    ds_nvect=0;
    Vector_3 k;

    for(l=0;l<nk;l++){ //all possible vectors
      for(i=0;i<3;i++){
	k[i]=(double)arr[MAXLEV*l+i];
      }
      k*=2*M_PI/anL;

      // determination of nvect
      int e;
      int m=0;
      while(k[m]==0. && m<3)m++;

      for(e=0;e<4;e++){ // averaging through 8 possible directions

	ds_kvect[ds_nvect++]=k;
	if(e!=0)k[m]=-k[m];

	m=(m+1)%3;
	if(k[m]==0.){
	  e+=2;
	  m=(m+1)%3;
	  if(k[m]==0.)break;
	}
	k[m]=-k[m];
      }
    }

  }

  printf("%d vectors\n",ds_nvect);
  if(ds_nvect==0)serror("!!!!!!???????\n");

  ds_e = new pfarr[ds_nvect];
  ds_p = new pfarr[ds_nvect];
  eold= new comple[ds_nvect];
  pold= new comple[ds_nvect];
  if(!ds_e || !ds_p || !eold || !pold)serror("DSF: MAE\n");

  for(i=0;i<ds_nvect;i++){
    eold[i][0]=eold[i][1]=pold[i][0]=pold[i][1]=0;
  }

  for(l=0;l<ds_nvect;l++){
    ds_e[l]= new double[2*ds_nf];
    ds_p[l]= new double[2*ds_nf];
    if(!ds_e[l] || !ds_p[l])serror("DSF: MAE\n");
  }

  dre= new comple[ds_nvect];
  drp= new comple[ds_nvect];
  if(!dre || !drp)serror("DSF: MAE\n");

  ds_count=0;
  ds_stpx=0;
  ds_stpi=0;
}




void Dstruc_step(){

  int j,m,l;
  int nn2=npart/2;
  Vector_3 r,k;


  for(l=0;l<ds_nvect;l++){ // assembling data for  all possible vectors
    k=ds_kvect[l];

    dre[l].x=dre[l].y=drp[l].x=drp[l].y=0.;

    for(j=0;j<npart;j++){ // averaging through nn2 particles
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];

      if(j<nn2){
	drp[l].x+=cos(r*k);
	drp[l].y+=sin(r*k);
      }
      else{
	dre[l].x+=cos(r*k);
	dre[l].y+=sin(r*k);
      }
    }
  }



  double d;
  int dumm;
  comple vr;

  while(ds_stpi<ds_nf && (double)ds_count>ds_stpx){

    d=DxScale(1.,stpx,dumm);

    for(l=0;l<ds_nvect;l++){
      for(m=0;m<2;m++)vr[m]=eold[l][m]*d+dre[l][m]*(1-d);
      ds_e[l][2*ds_stpi]  =  vr[0]/nn2;
      ds_e[l][2*ds_stpi+1]=  vr[1]/nn2;

      for(m=0;m<2;m++)vr[m]=pold[l][m]*d+drp[l][m]*(1-d);
      ds_p[l][2*ds_stpi]  =  vr[0]/nn2;
      ds_p[l][2*ds_stpi+1]=  vr[1]/nn2;
    }

    ds_stpi++;
    ds_stpx=((double)(nsteps-1))*ds_stpi/(ds_nf-1);

  }

  for(l=0;l<ds_nvect;l++){
    for(m=0;m<2;m++){
      pold[l][m]=drp[l][m];
      eold[l][m]=dre[l][m];
    }
  }

  ds_count++;
}



void Dstruc_process(char *dname){
  char str[256]="dynstruc.dat";

  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-dd.dat");
  int i,l;

  FILE *ff=fopen("tmp.dat","wt");
  for(i=0;i<ds_nf;i++){
    double ar=0,ai=0;
    for(l=0;l<ds_nvect;l++){
      ar+=ds_e[l][2*i];
      ai+=ds_e[l][2*i+1];
    }
    ar/=ds_nvect;
    ai/=ds_nvect;
    fprintf(ff,"%e %e %e\n",i*dt*Wp,ar,ai);
  }
  fclose(ff);


  FILE *f1=Err_fopen(str,"wt");
  printf("Writing dynamical structure factor...\n");
  double f=ds_f0;
  double t0=(nsteps*dt)/ds_nf;

  if(!time_aver){
    double t;
    fprintf(f1,"#1-t 2-Re(See) 3-Im(See) 4-Re(Spp) 5-Im(Spp) 6-Re(Spe) "
               "7-Im(Spe) 8-Re(Sep) 9-Im(Sep) 10-Re(S) 11-Im(S)\n");
    double aeer,aeei,appr,appi,aepr,aepi,aper,apei,Sr,Si;
    for(i=0;i<ds_stpi;i++){
     t= i*t0;
     aeer=aeei=appr=appi=aepr=aepi=aper=apei=0.;
     for(l=0;l<ds_nvect;l++){
     		aeer+=ds_e[l][0]*ds_e[l][2*i]+ds_e[l][1]*ds_e[l][2*i+1];
       aeei+=ds_e[l][2*i]*ds_e[l][1]-ds_e[l][2*i+1]*ds_e[l][0];

       appr+=ds_p[l][0]*ds_p[l][2*i]+ds_p[l][1]*ds_p[l][2*i+1];
       appi+=ds_p[l][2*i]*ds_p[l][1]-ds_p[l][2*i+1]*ds_p[l][0];

       aper+=ds_p[l][0]*ds_e[l][2*i]+ds_p[l][1]*ds_e[l][2*i+1];
       apei+=ds_e[l][2*i]*ds_p[l][1]-ds_e[l][2*i+1]*ds_p[l][0];

       aepr+=ds_e[l][0]*ds_p[l][2*i]+ds_e[l][1]*ds_p[l][2*i+1];
       aepi+=ds_p[l][2*i]*ds_e[l][1]-ds_p[l][2*i+1]*ds_e[l][0];
     }
     appr/=ds_nvect;
     aeer/=ds_nvect;
     appi/=ds_nvect;
     aeei/=ds_nvect;
     aepr/=ds_nvect;
     aepi/=ds_nvect;
     aper/=ds_nvect;
     apei/=ds_nvect;

     Sr=appr+aeer-aper-aepr;
     Si=appi+aeei-apei-aepi;



     fprintf(f1,"%e %e %e %e %e %e %e %e %e %e %e\n",t*Wp,aeer,aeei,
                                           appr,appi,aper,apei,aepr,aepi,
                                           Sr,Si);

     f+=ds_fc;
    }
    fclose(f1);
    return;
  }



  if(ds_dir){
    double *tmp;
    for(i=0;i<ds_nvect;i++){
      tmp=four_direct(ds_e[i],ds_nf,1,t0*ds_f0*Wp,t0*ds_f1*Wp,t0*ds_fc*Wp);
      delete ds_e[i];
      ds_e[i]=tmp;
      tmp=four_direct(ds_p[i],ds_nf,1,t0*ds_f0*Wp,t0*ds_f1*Wp,t0*ds_fc*Wp);
      delete ds_p[i];
      ds_p[i]=tmp;
    }
    ds_nf=2*((int)((ds_f1-ds_f0)/ds_fc)-1);
  }
  else{
    for(i=0;i<ds_nvect;i++){
      four(ds_e[i],ds_nf,1);
      //fold_it(ds_e[i],ds_nf);
      four(ds_p[i],ds_nf,1);
      //fold_it(ds_p[i],ds_nf);
    }
    ds_f0=0;
    ds_fc=2*M_PI/(nsteps*dt*Wp);
  }

  f=ds_f0;
  fprintf(f1,"#1-Wpe 2-Spp 3-See 4-Re(Sep) 5-Im(Sep) 6-Re(S) 7-Im(S)\n");
  double t,aee,app,aepr,aepi,S,Sr;
  for(i=0;i<=ds_nf/2;i++){
    t= i*dt;
    aee=app=aepr=aepi=0;
    for(l=0;l<ds_nvect;l++){
      app+= ds_p[l][2*i]*ds_p[l][2*i]+ds_p[l][2*i+1]*ds_p[l][2*i+1];
      aee+= ds_e[l][2*i]*ds_e[l][2*i]+ds_e[l][2*i+1]*ds_e[l][2*i+1];
      aepr+=ds_p[l][2*i]*ds_e[l][2*i]+ds_p[l][2*i+1]*ds_e[l][2*i+1];
      aepi+=ds_p[l][2*i+1]*ds_e[l][2*i]-ds_p[l][2*i]*ds_e[l][2*i+1];
    }

    app/=ds_nvect;
    aee/=ds_nvect;
    aepr/=ds_nvect;
    aepi/=ds_nvect;

    Sr=aee+app-2*aepr;
    S=Sr; //sqrt(Sr*Sr+4*aepi*aepi);
    fprintf(f1,"%e %e %e %e %e %e %e\n",
	    f,app, aee,aepr,aepi,S,-2*aepi);
    f+=ds_fc;
  }
  fclose(f1);
}




int nd_n, nd_type, nd_nvect;
Vector_3 *nd_kvect;

typedef Correlation *CorrP;


# define NCORR 5  // 0-i, 1-e, 2-ie, 3-ei, 4-sum (for check only)
CorrP nd_r[NCORR];


void NDstruc_init(){
  char tmpstr[250];
  int l,i;

  if(!(tType&COORD))serror("Can't calculate dynamic structure factor,\n"
			   " trajectory file does not contain Coords data.\n");
  unsigned int arr[MAXLEV*MAXSPLITS];

  Set_stop(1);

  Read_param("DSF calculation type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))nd_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))nd_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))nd_type=VC_STRAIT;
  else serror("New DSF: unknown calculation type.\n");

  if(nd_type==VC_FFT){
    printf("New DSF: Fourier degree: %d",(int)log2(nsteps)+1);
    nd_n=1<<((int)log2(nsteps)+1);
    printf(", %d points\n",vc_n);
  }
  else nd_n=nsteps;

  Set_stop(0);

  double kval;
  if(Read_param("New DSF Free K: %lf %d", &kval, &nd_nvect)==2){
    printf("Warning: using not allowed k-directions! \n");
    nd_kvect= new Vector_3[nd_nvect];
    if(!nd_kvect)serror("New DSF: MAE\n");
    for(i=0;i<nd_nvect;i++){
      nd_kvect[i]=Vector_3(random2(),random2(),random2());
      nd_kvect[i].normalize();
      nd_kvect[i]*=kval*2*M_PI/anL;
    }
  }
  else{

    Set_stop(1);
    int ksq;
    Read_param("New DSF K square: %d",&ksq);
    //Read_param("Frequency interval, steps: %f, %f, %d",&ds_f0,&ds_f1,&ds_nf);

    int nk;
    nk=splitsq(ksq,3,arr);
    printf("%d components\n",nk);



    nd_kvect= new Vector_3[8*nk];
    if(!nd_kvect)serror("DSF: MAE\n");

    nd_nvect=0;
    Vector_3 k;

    for(l=0;l<nk;l++){ //all possible vectors
      for(i=0;i<3;i++){
	      k[i]=(double)arr[MAXLEV*l+i];
      }
      k*=2*M_PI/anL;

      // determination of nvect
      int e;
      int m=0;
      while(k[m]==0. && m<3)m++;

      for(e=0;e<4;e++){ // averaging through 8 possible directions

       nd_kvect[nd_nvect++]=k;
       if(e!=0)k[m]=-k[m];

       m=(m+1)%3;
       if(k[m]==0.){
         e+=2;
         m=(m+1)%3;
         if(k[m]==0.)break;
       }
       k[m]=-k[m];
      }
    }

  }

  printf("%d vectors\n",nd_nvect);
  if(nd_nvect==0)serror("!!!!!!???????\n");

  int m;
  for(m=0;m<NCORR;m++){
   nd_r[m] = new Correlation[nd_nvect];
   if(!nd_r[m])serror("New DSF: MAE\n");

   for(i=0;i<nd_nvect;i++){
    nd_r[m][i].init(nd_n, nsteps*dt*Wp,CORR_COMPLEX);
    if(m<2 || m>3)nd_r[m][i].insert_begin(nsteps); // only for ee and ii
   }
  }

}




void NDstruc_step(){

  int j,l;

  Vector_3 r,k;
  double sr,si;

  for(l=0;l<nd_nvect;l++){ // assembling data for  all possible vectors
    k=nd_kvect[l];

    sr=si=0.;
    for(j=0;j<nion;j++){ // averaging through nn2 particles
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];

      sr+=cos(r*k);
      si+=sin(r*k);
    }
    nd_r[0][l].insert_next(sr,si);
    double a=sr, b=si;

    sr=si=0.;
    for(j=nion;j<npart;j++){ // averaging through nn2 particles
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];

      sr+=cos(r*k);
      si+=sin(r*k);
    }
    nd_r[1][l].insert_next(sr,si);

    nd_r[4][l].insert_next(a-sr,b-si);

  }

}



void NDstruc_process(char *dname){
  char str[256]="dynstruc.dat";

  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-dy");
  int i;

  for(i=0;i<nd_nvect;i++){
   nd_r[0][i].insert_end();
   nd_r[1][i].insert_end();
  }


  printf("Writing NEW dynamical structure factor...\n");

  if(nd_type==VC_DIRECT){
   //strcat(filename,"d.dat");
   strcat(str,".dat");
   printf("New DSF: Direct correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     printf("Vector #%d\n",i+1);
     CrossCorrDirect(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].direct();
   }
  }
  else if(nd_type==VC_STRAIT){
   //strcat(filename,"s.dat");
   strcat(str,".dat");
   printf("New DSF: Strait correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     printf("Vector #%d\n",i+1);
     CrossCorrStrait(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].strait();
   }
  }
  else{ //FFT
   //strcat(filename,"f.dat");
   strcat(str,".dat");
   printf("New DSF: FFT correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     printf("Vector #%d\n",i+1);
     CrossCorrFFT(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].calculate();
   }
  }

  // averaging through vectors
  int n=nd_r[0][0].n;
  int j,m;
  for(j=0;j<n;j++){
   for(m=0;m<NCORR;m++){
    for(i=1;i<nd_nvect;i++){
     nd_r[m][0].arr[j]+=nd_r[m][i].arr[j]; // real part
     nd_r[m][0].arri[j]+=nd_r[m][i].arri[j]; // imaginary part
     if(nd_type==VC_FFT){
      nd_r[m][0].arr[n+j]+=nd_r[m][i].arr[n+j]; // fft^2
      // only for cross-correlations
      if(m>=2)nd_r[m][0].arri[n+j]+=nd_r[m][i].arri[n+j]; // imaginary of fft1 mal fft2*
     }
    }
    nd_r[m][0].arr[j]/=n;
    nd_r[m][0].arri[j]/=n;

    if(nd_type==VC_FFT){
     nd_r[m][0].arr[n+j]/=n;
     if(m>=2)nd_r[m][0].arri[n+j]/=n;
    }
   }
  }

  FILE *f1=Err_fopen(str,"wt");

  //for(i=0;i<NCORR;i++){
  // nd_r[i][0].normalize();
  //}

  double Sr, Si;
  double tmax=nd_r[0][0].maxtime;
  double dt=0.1,t;

  fprintf(f1,"#1-t 2-Re(Sii) 3-Im(Sii) 4-Re(See) 5-Im(See) 6-Re(Sie) "
               "7-Im(Sie) 8-Re(Sei) 9-Im(Sei) 10-Re(S) 11-Im(S)\n");

  for(t=0.;t<tmax;t+=dt){
   Sr=nd_r[0][0].r(t)+nd_r[1][0].r(t)-nd_r[2][0].r(t)-nd_r[3][0].r(t);
   Si=nd_r[0][0].i(t)+nd_r[1][0].i(t)-nd_r[2][0].i(t)-nd_r[3][0].i(t);
   fprintf(f1,"%f %e %e %e %e %e %e %e %e %e %e %e %e\n",t,nd_r[0][0].r(t),nd_r[0][0].i(t),
                                                     nd_r[1][0].r(t),nd_r[1][0].i(t),
                                                     nd_r[2][0].r(t),nd_r[2][0].i(t),
                                                     nd_r[3][0].r(t),nd_r[3][0].i(t),
                                                     Sr,Si,
                                                     nd_r[4][0].r(t),nd_r[4][0].i(t));
  }

  fclose(f1);
  if(nd_type==VC_FFT){
   str[strlen(str)-4]=0;
   strcat(str,"f.dat");
   f1=Err_fopen(str,"wt");

   fprintf(f1,"#1-w 2-Fii^2 3-Fee^2 4-Re(Fie) 5-Im(Fie) 6-Re(Fei) "
               "7-Im(Fei) 8-Re(F^2) 9-Im(F^2)\n");

   double w,wmax=10.,dw=0.1;
   for(w=0.;w<wmax;w+=dw){
    Sr=nd_r[0][0].four_sq(w)+nd_r[1][0].four_sq(w)-nd_r[2][0].four_sq(w)-nd_r[3][0].four_sq(w);
    Si=-nd_r[2][0].four_sqi(w)-nd_r[3][0].four_sqi(w);
    fprintf(f1,"%f %e %e %e %e %e %e %e %e %e %e\n",w,nd_r[0][0].four_sq(w),
                                                nd_r[1][0].four_sq(w),
                                                nd_r[2][0].four_sq(w),
                                                nd_r[2][0].four_sqi(w),
                                                nd_r[3][0].four_sq(w),
                                                nd_r[3][0].four_sqi(w),
                                                Sr,Si,
                                                nd_r[4][0].four_sq(w),nd_r[4][0].four_sqi(w));
   }
   fclose(f1);
  }

  for(i=0;i<NCORR;i++){
   delete[] nd_r[i];
  }

}





extern Distribution DRRee,DRRep,DRRpp;


void Pair_init(){
  if(!(tType&COORD))serror("Can't calculate pair distribution function,\n"
			   " trajectory file does not contain Coords data.\n");

  Set_stop(0);
  double rr_r0, rr_r1;
  if(Read_param("r-r range: %lf, %lf",&rr_r0,&rr_r1)!=2){
    rr_r0=0.;
    rr_r1=-1.;
  }

  DRRee.init(rr_r0,(rr_r1>0 ? rr_r1 : anL/2),400);
  DRRep.init(rr_r0,(rr_r1>0 ? rr_r1 : anL/2),400); 
  if(non_symm)DRRpp.init(0,(rr_r1>0 ? rr_r1 : anL/2),400);
}



void Pair_step(){
  int i,j;
  double R;
  Vector_3 ri,r;  // very ugly!

  int ni=npart/2; // only for charge-symmetrical case!

  int type=0; //IONION

  for(i=0;i<npart;i++){
    if(i>=ni)type|=0x1;// setting electron-? interaction type
    type&=0x1;         // setting ?-ion interaction type

    if(!xx_used)ri=r_cell(anx[i]);
    else ri=anx[i];

    for(j=i+1;j<npart;j++){
 
      if(!xx_used)r=r_cell(anx[j]);
      else r=anx[j];
      
      r=r-ri;
     
      int k;
      for(k=0;k<3;k++){
	if(r[k]>anL/2)r[k]-=anL;
	if(r[k]<-anL/2)r[k]+=anL;
      }

      R=r.norm();
      
      if(type==ELCELC || type ==IONION){      
	if(non_symm && type==IONION)DRRpp.point(R,1.);
	else DRRee.point(R,1.);
      }
      else DRRep.point(R,1.);
    }
  }
}



void WriteDRR(char *file){
FILE *f1=Err_fopen(file,"wt");
 int i;
 double a,b,c,x,dx,aa,bb,cc;
 double l=(DRRee.x2-DRRee.x1); 
 dx=(DRRee.x2-DRRee.x1)/399;

 //xtest=DRRee.x2-3*dx;
 //double na=DRRee(xtest);
 //double nb=DRRep(xtest);
 //double nrma=-4*M_PI*(RDebye*RDebye)/(DRRee.norm-na*l*l*l/(3*xtest*xtest));
 //double nrmb=+4*M_PI*(RDebye*RDebye)/(DRRep.norm-nb*l*l*l/(3*xtest*xtest));
 //double intg=1./(xtest*xtest);

 double na=0.;
 if(DRRee.count)na=l*l*399/3/DRRee.count;
 double nb=0.;
 if(DRRep.count)nb=l*l*399/3/DRRep.count;
 double nc=0.;
 if(non_symm && DRRpp.count)nc=l*l*399/3/DRRpp.count;


 for(x=DRRee.x1,i=0;i<400;x+=dx,i++){
   a=DRRee(x);
   b=DRRep(x);
   if(fabs(x)>1e-10){
    aa=a*na/(x*x);
    bb=b*nb/(x*x);
   }
   else aa=bb=0.;

   if(!non_symm){
     fprintf(f1,"%f %f %f %f %f \n",x,
	   //a,b,(a-na*x*x*intg)*nrma,(b-nb*intg*x*x)*nrmb, 4*M_PI*x*exp(-x/RDebye));

	     a*na,b*nb,aa,bb
	     /*, 4*M_PI*x*exp(-x/anRDebye)*/);
   }
   else{
     c=DRRpp(x);
     if(fabs(x)>1e-10)cc=c*nc/(x*x);
     else cc=0.;
   
     fprintf(f1,"%f %f %f %f %f %f %f \n",x,
     	     a*na,c*nc,b*nb,aa,cc,bb
	     /*, 4*M_PI*x*exp(-x/anRDebye)*/);
   }

 }
 fclose(f1);
}


void Pair_process(char *dname){
  char str[256]="statstruc.dat";
  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-rr.dat");
  WriteDRR(str);
}


Distribution DVe[4],DVp[4];


void V_init(){
  if(!(tType&VEL))serror("Can't calculate velocity distribution function,\n"
			 "trajectory file does not contain Velocity data.\n");
 
  Set_stop(0);
  double ve0=0.,ve1=2.;
  double vp0=0., vp1=2./sqrt(animass);
  int gvp=400, gve=400;
  Read_param("Ve distribution range: %lf, %lf",&ve0,&ve1);
  Read_param("Ve distribution grid: %d",&gve);
  Read_param("Vi distribution range: %lf, %lf",&vp0,&vp1);
  Read_param("Vi distribution grid: %d",&gvp);


  int i;
  for(i=0;i<3;i++){
    DVe[i].init(ve0,ve1,gve);
    DVp[i].init(vp0,vp1,gvp);
  }
  DVe[3].init(0.,ve1,gve);
  DVp[3].init(0.,vp1,gvp);
}



void V_step(){
  int i,j;
  for(i=0;i<nion;i++){
    for(j=0;j<3;j++){
      DVp[j].point(anv[i][j],1.);
    }
    DVp[3].point(anv[i].norm(),1.);
  }
  for(i=nion;i<npart;i++){
    for(j=0;j<3;j++){
      DVe[j].point(anv[i][j],1.);
    }
    DVe[3].point(anv[i].norm(),1.);
  }
}




void V_process(char *dname){
  char str[256]="veldistr.dat";



  int i,j;
  double x,dx;
  for(i=0;i<4;i++){
    DVp[i].normalize((double)DVp[i].count/DVp[i].allcount);
    DVe[i].normalize((double)DVe[i].count/DVe[i].allcount);
  }




  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-vi.dat");

  FILE *f1=Err_fopen(str,"wt");
  dx=(DVp[0].x2-DVp[0].x1)/(DVp[0].n-1);
  for(x=DVp[0].x1,i=0;i<DVp[0].n;i++,x+=dx){
    fprintf(f1,"%f  ",x);
    for(j=0;j<3;j++)fprintf(f1,"%f  ",DVp[j](x));
    if(x<0)fprintf(f1,"%f  ",0.);
    else fprintf(f1,"%f  ",DVp[3](x));
    fprintf(f1,"\n");
  }
  fclose(f1);


  strcpy(str,out_dir);
  strcat(str,dname);
  strcat(str,"-ve.dat");
  f1=Err_fopen(str,"wt");
  dx=(DVe[0].x2-DVe[0].x1)/(DVe[0].n-1);
  for(x=DVe[0].x1,i=0;i<DVe[0].n;i++,x+=dx){
    fprintf(f1,"%f  ",x);
    for(j=0;j<3;j++)fprintf(f1,"%f  ",DVe[j](x));
    if(x<0)fprintf(f1,"%f  ",0.);
    else fprintf(f1,"%f  ",DVe[3](x));
    fprintf(f1,"\n");
  }
  fclose(f1);

}






int cur_cur=0, refine=0,pairdistr =0, dynstruc=0,vdistr=0;
int vcorr=0;
int dmoment=0;
int ndstruct=0;

char dname[100];
char gdname[100];
double max_time=-1., start_time=0.;
int trjstep=1;


void init_all(){
  if(refine) Refine_init(dname);
  if(micro_field) Microfield_init();
  if(cur_cur)Cur_init();
  if(statstruc)Sstruc_init();
  if(pairdistr)Pair_init();
  if(dynstruc)Dstruc_init();
  if(vdistr)V_init();
  if(vcorr)Vcorr_init();
  if(dmoment)DM_init();
  if(ndstruct)NDstruc_init();
}




void read_cfg(){
  char tmpstr[256];
  //Read_param("Trajectory file: %s",datafile);

  Set_stop(0);

  if(Read_param("Microfield analysis: %s",tmpstr)){
    if(tmpstr[0]=='y')micro_field=1;
  }

  if(Read_param("Current-current corr*: %s",tmpstr)){
   if(tmpstr[0]=='y')cur_cur=1;
  }

  if(Read_param("Refine: %s",tmpstr)){
    if(tmpstr[0]=='y')refine=1;
  }

  if(Read_param("Structure factor: %s",tmpstr)){
    if(tmpstr[0]=='y')statstruc=1;
  }


  if(Read_param("Pair distribution function: %s",tmpstr)){
    if(tmpstr[0]=='y')pairdistr=1;
  }

  if(Read_param("Dynamic structure factor: %s",tmpstr)){
    if(tmpstr[0]=='y')dynstruc=1;
  }

  if(Read_param("Velocity distribution function: %s",tmpstr)){
    if(tmpstr[0]=='y')vdistr=1;
  }

  time_aver=1;
  if(Read_param("Time averaging: %s",tmpstr)){
    if(tmpstr[0]=='n')time_aver=0;
  }

  vcorr=1;
  if(Read_param("Velocity autocorrelation: %s",tmpstr)){
    if(tmpstr[0]=='n')vcorr=0;
  }

  dmoment=1;
  if(Read_param("Dipole moment: %s",tmpstr)){
    if(tmpstr[0]=='n')dmoment=0;
  }

  ndstruct=1;
  if(Read_param("New DSF: %s",tmpstr)){
    if(tmpstr[0]=='n')ndstruct=0;
  }

# ifdef UNIX

  if(Read_param("Output directory: %s",tmpstr)){
    strcpy(out_dir,tmpstr);
    if(out_dir[strlen(out_dir)]!='/')strcat(out_dir,"/");

    if(mkdir(out_dir,S_IRWXU|S_IRGRP|S_IROTH)){
      if(errno!=EEXIST)serror("Can not create directory: %s.\n%s\n",
			     out_dir,strerror(errno));
    }

  }
# endif

  if(!Read_param("Stop  time: %lf",&max_time)){
    if(Read_param("Stop  step: %lf",&max_time))
      max_time*=dt*Wp;
  }
  if(!Read_param("Start time: %lf",&start_time)){
    if(Read_param("Start step: %lf",&start_time))
      start_time*=dt*Wp;
  }

  max_time/=Wp;
  start_time/=Wp;
  Read_param("Trajectory step: %d",&trjstep);
  if(max_time<0)max_time=dt*nsteps;
  if(max_time<start_time)serror("Start and/or stop time are not consistent with trajectory length.\n");

  nsteps=(long)((max_time-start_time)/dt);
  nsteps/=trjstep;
  dt*=trjstep;

  Read_param("Set data name: %s",dname);

  strcpy(gdname,dname);
  strcat(gdname,"g.an");

}

void one_step(){

 if(BadStatus)printf("Reading procedure returned 'bad status'...\n");

 if(cur_step*dt<start_time){
   return;
 }
 if(max_time>0 && cur_step*dt>max_time){
  cur_step--;
  return;
 }
 if(refine)Refine_step();
 if(cur_cur)Cur_step();
 if(micro_field)Microfield_step();
 if(statstruc)Sstruc_step();
 if(pairdistr)Pair_step();
 if(dynstruc)Dstruc_step();
 if(vdistr)V_step();
 if(vcorr)Vcorr_step();
 if(dmoment)DM_step();
 if(ndstruct)NDstruc_step();
 return;
}

int Init_analyse(char *cfgfile,int nn,double L, double Gamma, double T,char *daname,double tstep, long nstp){

  if(!cfgfile || strstr(cfgfile,"default"))cfgfile="analyse.cfg";
  strcpy(dname,daname);
  anL=L;
  anGamma=Gamma;
  anT=T;
  npart=nn;
  nion=nn/2;
  nelec=nn/2;
  anRDebye=1./sqrt(4*M_PI*npart/(anL*anL*anL));
  dt=tstep;
  tType=VEL|COORD;
  Open_param_file(cfgfile);
  read_cfg();
  init_all();
  Close_param_file();
  nsteps=nstp;
  cur_step=0;
  return 0;
}




void Step_analyse(Vector_3 *xx, Vector_3 *v){
  xx_used=1;
  anx=xx;
  anv=v;

  one_step();
  //nsteps++;
// cur_step++;
}

void Process_analyse(){
  if(refine)Refine_process();
  if(cur_cur)Cur_process(dname);
  if(micro_field)Microfield_process(dname);
  if(statstruc)Sstruc_process(dname);
  if(pairdistr)Pair_process(dname);
  if(dynstruc)Dstruc_process(dname);
  if(vdistr)V_process(dname);
  if(vcorr)Vcorr_process(dname);
  if(dmoment)DM_process(dname);
  if(ndstruct)NDstruc_process(dname);
}














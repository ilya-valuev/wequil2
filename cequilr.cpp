#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <fcntl.h>

# ifdef UNIX

#include <sys/stat.h> // define mkdir only
#include <sys/types.h>
#include <unistd.h> // defines execl only
#include <errno.h> // defines errno only

//# define va_start 234

# else

# include <io.h>
# include <direct.h>
# include <errno.h>

# endif

#include "vector_3.h"
#include "common.h"
#include "rparam.h"
#include "statist.h"
#include "analyse.h"
#include "pcycle.h"
#include "interact.h"



# include "plasma.h"
# include "plasmar.h"


double Efm, Tr;
Vector_3 RandomForce(Vector_3 & v){
 static Vector_3 f;
 for(int k=0;k<3;k++){
   f[k]=Efm*random1()-Tr*v[k];
   //Efm*2*(0.5-(double)rrandom()/RAND_MAX)-Tr*v[k];
 }
 return f;
}

double getEmax(Plasma& Gas){
  int type=ELCION;
  double dEpot, dEquant;
  Gas.potential(type, 1e-12,dEpot,dEquant);
  return fabs(dEpot);
}



class mc_simm{
public:
  int momentum_depend;
  int auto_adjust;
  Vector_3 *newc;
  double E1;

  int n;
  double dc[2];
  long accept[2];
  long aold[2];
  long nav;
  long nstp;
  double T;
  double Epotent;
  int one_particle;
  int diff_temp;
  mc_simm(Plasma &Gas, double d1,double d2, int is_one=0, int is_diff=0){
    momentum_depend=0;
    auto_adjust=1;
    one_particle=is_one;
    diff_temp=is_diff;
    n=Gas.n;

    dc[1]=d1;
    dc[0]=d2;
    accept[0]=accept[1]=0;
    aold[0]=aold[1]=0;
    nstp=0;
    nav=50;

    if(one_particle){
      Gas.init_matrix();
      newc= new Vector_3;
    }
    else newc= new Vector_3[n];
    if(!newc)serror("mc_simm: Memory alocation error.\n");

    T=Gas.getT();
    Epotent=Gas.Epotent;
    if(!diff_temp)E1=1.5*T*n+Epotent;
    else E1=1.5*Gas.Ti*Gas.ni+1.5*Gas.Te*Gas.ne+Epotent;
  }
  ~mc_simm(){
    delete newc;
  }


  // ratio of accepted/all steps
  double get_ratio(int sttype=1){
    if(!nstp) return 0;
    else return ((double)(accept[sttype]+aold[sttype]))/((double)nstp);
  }

  void initial(Plasma &Gas){
    T=Gas.getT();
    Gas.interaction();
    Epotent=Gas.Epotent;
    if(!diff_temp)E1=1.5*T*n+Epotent;
    else E1=1.5*Gas.Ti*Gas.ni+1.5*Gas.Te*Gas.ne+Epotent;
  }

  int step(Plasma &Gas);
}  *MC;

int mc_simm::step(Plasma &Gas){
  int i,j;


  Vector_3 *tmp=newc;

  static int sttype=1;

  int m;
  if(one_particle){
    m=(int)(Gas.n*(0.5+random1()));
    if(m>=Gas.n)m=Gas.n-1;
    if(sttype){
      newc[0]=Gas.x[m];
      for(j=0;j<3;j++){
      	Gas.x[m][j]+=dc[sttype]*random1();
      }
    }
    else{
      newc[0]=Gas.v[m];
      for(j=0;j<3;j++){
       	Gas.x[m][j]+=dc[sttype]*random1();
      }
    }
  }
  else{
    m=-1;
    if(sttype)tmp=Gas.x;
    else tmp=Gas.v;
    for(i=0;i<Gas.n;i++){
      for(j=0;j<3;j++){
        	newc[i][j]=tmp[i][j]+dc[sttype]*random1();
      }
    }
    if(sttype)Gas.x=newc;
    else Gas.v=newc;
  }


  int mult=1;
  double Epot_old=Gas.Epotent;
  if(momentum_depend){
    Gas.interaction(m);
    T=Gas.getT();
    Epotent=Gas.Epotent;
    mult=2;
  }
  else{
    if(sttype){
      Gas.interaction(m,1);
      Epotent=Gas.Epotent;
    }
    else T=Gas.getT();
  }

  double E2;
  if(!diff_temp)E2=1.5*T*Gas.n+Epotent;
  else E2=1.5*Gas.Ti*Gas.ni+1.5*Gas.Te*Gas.ne+Epotent;

  double Temp=1.;
  if(diff_temp && m>=0){
    if(m<Gas.ni)Temp=Gas.ini_Ti;
    else Temp=Gas.ini_Te;
  }


  int ac =0;
  if(E2<=E1)ac=1;
  else{
     if(exp((E1-E2)/Temp)>0.5+random1())ac=1;
  }
  (accept[sttype])+=ac;
  nstp++;

  if(nstp%nav==0){
    aold[0]=accept[0];
    aold[1]=accept[1];
    accept[1]=accept[0]=0;
    nstp=nav;
  }


  double ratio=mult*get_ratio(sttype);

  if(auto_adjust){
   double k;
   double kmax=1.1,kmax_inv=1./kmax;
   if((nstp/2)%5==0){
     if(ratio<0.5*kmax_inv)k=kmax_inv;
     else if(ratio>0.5*kmax)k=kmax;
     else k=ratio/0.5;

     dc[sttype]*=k;
     if(sttype && dc[sttype]>Gas.L/2)dc[sttype]=Gas.L/2;
   }
   //printf("%d, %f, %f\n",(one_particle ? m : sttype),ratio,dc[sttype]);
  }

  if(ac){
    if(!one_particle)newc=tmp;
    E1=E2;
  }
  else{ // restoring previous
    if(one_particle){

      for(i=0;i<Gas.n;i++){
	      (*Gas.umatr)(i,m)=(*Gas.umatr)(i,i);
	       //printf("%i -> (%d, %d) [%f]\n",i,i,m,fucktmp[i]);
      }
      if(sttype){
	      Gas.x[m]=newc[0];
       Gas.rcell(m);
      }
      else Gas.v[m]=newc[0];

      //Gas.Epotent=Epot_old;
      //Gas.interaction(m);
      // printf("%f\n",Gas.Epotent-Epot_old);

    }
    else{
      if(sttype)Gas.x=tmp;
      else Gas.v=tmp;
    }

    Gas.Epotent=Epot_old;
    //Gas.Epotent=E1-1.5*T*Gas.n;
  }

  if(momentum_depend)sttype=1-sttype;

  /*
  if(nstp%4==0){
    //Gas.is_matr=0;
    Gas.interaction(m,1); //to check!
    //Gas.is_matr=1;
  }*/

  return ac;
}


int mc_equil=0, mc_calc=0;

char logfile[200];
double RFstatus=0.;
int lwrites=0;

void WriteToLog(double t,Plasma &Gas,Statistics *local){
  FILE *f1;
  f1=fopen(logfile,"at");
  if(!f1){
    printf("Can't write to log!\n");
    return;
  }
  if(!lwrites)// header
    fprintf(f1,"#1-t 2-Ekin 3-Epot 4-Etotal 5-Ecoul 6-Equant 7-mode\n");
  double Ekin,Epot,Efull,Ecoul,Equant, Ekini, Ekine;

  Ekin=1.5*local[0].av();
  if(Gas.stable_ions)Ekin*=((double)Gas.ne/Gas.n);
  Epot=local[2].av()/Gas.n;
  Efull=Ekin+Epot;
  Ecoul=local[1].av()/Gas.n;
  Equant=local[3].av()/Gas.n;

  Ekini=1.5*local[5].av();
  Ekine=1.5*local[6].av();

  fprintf(f1,"%f %f %f %f %f %f %f",
	  t,Ekin,Epot,Efull,Ecoul,Equant,RFstatus);
  if(mc_equil==1)fprintf(f1," %f %f\n",
		      (double)((double)(MC->accept[0]+MC->aold[0])/MC->nstp),
		      (double)((double)(MC->accept[1]+MC->aold[1])/MC->nstp));
  else fprintf(f1," %f %f\n",Ekini,Ekine);
  lwrites++;
  fclose(f1);
}

PlasmaRecord Trajectory;

extern int rel_step;
double Estab;
int catch_average=0;


int lim_sw=0;
int limrescale_switch(int sw){
  int t=lim_sw;
  lim_sw=sw;
  return t;
}

double limTe=-1, limTi=-1, limT=-1;
int limstpe=0,limstpi=0,limstp, limspec=0;;
long gcount=0;

double fixT=-1;
double incT0=-1,  incdT=-1;
long incStp=0, incEnd=0;
double lim_t=0;

double MoveIt(double &t,long wr_int,long nsteps,Plasma& Gas,Statistics** stats,Statistics *local,int nstat){


 long i;
 int j;

 Statistics *temp= new Statistics[nstat];
 if(!temp)fatal_error("MoveIt: memory allocation error\n");

 double E0=Gas.Epotent;

 for(j=0;j<nstat;j++){
  temp[j]=*(stats[j]);
  temp[j].clear();
  local[j].clear();
 }
 int tmp_wrdistr=write_distr;

 for(i=0;i<nsteps;i++){

  if(i%wr_int ==0)write_distr=tmp_wrdistr; // activating RR-distr write
  else write_distr=0;

  int ac=1;
  if(mc_equil==1){
    ac=MC->step(Gas);
    ac=1;
  }
  else{
    Gas.interaction();
    Gas.stepmove();
  }

  for(j=0;j<nstat;j++){
    if(ac || i==0)temp[j].next();
    else temp[j].prev();
  }

  t+=Gas.Wpe*Gas.dt;

  if(i==0 && RFstatus!=0.){ // initial state write at every equillibration proc
    double tmp=RFstatus;
    RFstatus=-0.5;
    WriteToLog(t,Gas,temp);
    RFstatus=tmp;
  }

  if((i+1)%wr_int==0 || i>=nsteps-1){ // now writing log
   WriteToLog(t,Gas,temp);
   for(j=0;j<nstat;j++){
     local[j]+=temp[j];
     temp[j].clear();
   }
  }

  if(catch_average){
    if((E0-Estab)*(Estab-Gas.Epotent)>=0){
      catch_average=0;
      break;
    }
    E0=Gas.Epotent;
  }

  if(Trajectory.valid){  // writing output
    //printf("stp\n");
    Trajectory.Step();
    //printf("done\n");
    int q=Trajectory.Query();
    while(q>=0){
      switch(q){
      case PR_IONCOORDS:
	Trajectory.Fill(Gas.x);
	printf("Wrote ion coords.\n");
	break;
      case PR_IONVEL:
	Trajectory.Fill(Gas.v);
	printf("Wrote ion velocities.\n");
	break;
      case PR_ELCCOORDS:
	Trajectory.Fill(Gas.x+Gas.ni);
	printf("Wrote electron coords.\n");
	break;
      case PR_ELCVEL:
	Trajectory.Fill(Gas.v+Gas.ni);
	printf("Wrote electron velocities.\n");
	break;
      case PR_FLOW:
	{
	  Vector_3 flow[2];
	  Flow(Gas.n,Gas.v,flow,flow+1);
	  Trajectory.Fill(flow);
	  //printf("wf\n");
	}
	break;
      case PR_UNKNOWN:
	msg_error("Unknown frame fill request!\n");
	break;

      }
      q=Trajectory.Query();
    }
  }

  // fixing temperature reached some value
  if(fixT>0 && Gas.T>fixT && !incEnd){
     // entering "shelf" regime
     limspec=0x03;
     limTe=Gas.Te;
     limTi=Gas.Ti;
     incEnd=gcount+incStp;
     lim_sw=1;
     FILE *f=fopen("shelf.d","at");
     fprintf(f,"%f %f\n",t,t-lim_t);
     fclose(f);
  }
  if(incEnd && gcount>=incEnd){
     // leaving "shelf" regime
     incEnd=0;
     lim_sw=0;
     // parameters for the next "shelf"
     fixT+=incdT;
     lim_t+=incStp*Gas.dt;
  }



  // velocity rescaling
  if(lim_sw){
    if(limspec&0x02 && Gas.Te> limTe){// maybe need to rescale electron velocities
      if(gcount%limstpe==0){
         printf("Rescaling electron velocity from %f to %f!\n",Gas.Te, limTe);
         Gas.evel_scale(limTe);
      }
    }

    if(limspec&0x01 && limstpi && Gas.Ti> limTi){// maybe need to rescale ion velocities
      if(gcount%limstpi==0){
         printf("Rescaling ion velocity from %f to %f!\n",Gas.Ti, limTi);
         Gas.ivel_scale(limTi);
      }
    }

    if(limspec&0x04 && limstp && Gas.T> limT){// maybe need to rescale all velocities
      if(gcount%limstp==0){
         printf("Rescaling all velocities from %f to %f!\n",Gas.T, limT);
         Gas.vel_scale(limT);
      }
    }
  }
  gcount++;


 }
 for(j=0;j<nstat;j++)(*(stats[j]))+=local[j];

 write_distr=tmp_wrdistr;

 delete [] temp;
 return nsteps*Gas.dt;
}


int wr_film=0;
int flm_count=0;

void WriteFilm(char *file,Plasma &Gas){
  char str[200];
  FILE *f1;
  int i;
  sprintf(str,"%s.gnu",file);

  f1=Err_fopen(str,"at");
  fprintf(f1,"plot '%s.flm%04de', '%s.flm%04dp'\n# pause 1\n",file,flm_count,file,flm_count);
  fclose(f1);
  
  sprintf(str,"%s.flm%04de",file,flm_count++);
  f1=Err_fopen(str,"wt");
  for(i=0;i<Gas.n/2;i++)fprintf(f1,"%f  %f\n",Gas.xx[i][0]/Gas.L,Gas.xx[i][1]/Gas.L);
  fclose(f1);
  str[strlen(str)-1]='p';
  f1=Err_fopen(str,"wt");
  for(i=Gas.n/2;i<Gas.n;i++)fprintf(f1,"%f  %f\n",Gas.xx[i][0]/Gas.L,Gas.xx[i][1]/Gas.L);
  fclose(f1);
}




double Gf=10.;

double eq_dt;
long eq_nsteps,chk_nsteps,rf_nsteps, tst_nsteps, sw_nsteps;
int scale_vel=0;

int soft_step=1;
int rel_step=0;

double rf_dtcoeff=1.;
double rf_dt0=0.0005;

double delta=0.5;

double in_cs=-1.;

long wr_int=0;

int in_distr=ZEROVEL;

char filmdir[256]="film/";


int auto_rf=0, rf_sw_off=0;

char sfile[256]="cequil.st";


enum regimes{ 
  ACCEL=0,
  ROUGH=1,
  SCALE=2,
  SWITCH_OFF=3,
  TEST=4 
};

enum rf_types{
  OFF=-1,
  GF=0,
  AUTO=1,
  FLUC=2,
  EXP=3

};


int cur_rf, cur_regime,cur_iter;
long cur_nsteps;
double cur_dt,cur_dc;


void set_regime(int regime,double delta,double dt,long nsteps){
  cur_regime=regime;
  cur_iter=0;
  cur_nsteps=nsteps;
  cur_dt=(rel_step ? dt*rf_dtcoeff : dt);
  cur_dc=delta;
  cur_rf=auto_rf;

  if(mc_equil)mc_equil=1;

  if(regime==ACCEL){
    if(cur_rf==FLUC)cur_rf=AUTO;
  }
  else if(regime==SWITCH_OFF){
    cur_rf=EXP;
  }
  else if(regime==TEST){
    cur_rf=OFF;
    if(mc_equil && !mc_calc)
      mc_equil=-1;
  }
}



int e_stab=0;
double stab_acc=0.1;
int no_test=0;




double Equillibrium(double &t,Plasma &Gas,Statistics **stats, int nstat){
  struct{
    char s[20];
    double a;
  } RFst[5]= { {"ACCELERATION",0.75},
	       {"ROUGH",1.},
	       {"SCALE",1.5},
	       {"SWITCH OFF", 0.5},
	       {"TEST", 0.25}};

  int j;

  Gas.r0=0.05;
  Gas.init_config(in_cs,in_distr);
  if(mc_equil)MC->initial(Gas);


  double  T_av, dc1,dc2,Estb=0;

  double T_val=Gas.getT();

  int sw_n=10;

  dc2=delta*sqrt(1./Gas.n);

  if(scale_vel)dc1=0.5*sqrt(1./Gas.n);
  else dc1=dc2;

  Statistics *statsl= new Statistics[nstat];
  if(!statsl)fatal_error("Equilibrium: memeory allocation error\n");

  for(j=0;j<nstat;j++){
    statsl[j]=*(stats[j]);
    stats[j]->clear();
  }
  if(write_distr){
    DRRee.clear();
    DRRep.clear();
    if(non_symm)DRRpp.clear();
  }

  int iter=0;
  double Gf1=Gf,Gf0=Gf;

  Gas.ext_force=RandomForce;
  set_regime(ACCEL,dc1,rf_dt0,rf_nsteps);

  do{
    Gas.dt=cur_dt;

    switch(cur_rf){
    case AUTO:
      Gf0=Gf1=1./sqrt(Gas.dt*sqrt((double)Gas.n));
      break;
    case FLUC:
      {
         double dEp=statsl[2].dev()/Gas.n;
         double dEt=statsl[4].dev()/Gas.n;

         if(dEt<1e-10)Gf1*=50.;
         else Gf1=Gf1*sqrt(dEp/dEt);
      }
      break;
    case EXP:
      Gf1=Gf0*(1.-((double)cur_iter)/sw_n);      // exp(-5./sw_n);
      if(Gf1<0.)Gf1=0.;
      break;
    }


    Tr=Gf1*Gf1/(24.);  // adjusting RF and friction
    Efm=Gf1/sqrt(Gas.dt);

    RFstatus=RFst[cur_regime].a;

    if(StopStatus(sfile,0)){
      StopStatus(sfile,-1);
      serror("Program interrupted!\n");
    }

    MoveIt(t,wr_int,cur_nsteps,Gas,stats,statsl,nstat);
    iter++;
    cur_iter++;

    T_av=statsl[0].av();

    if(soft_step && cur_regime<SCALE)Gas.dt=cur_dt*sqrt(1./T_av);


    printf("%s: Average values for %ld steps with RF:\n"
	   "Epot=%f, Ecoul=%f, Temp=%f\n",RFst[cur_regime].s,cur_nsteps,
            statsl[2].av()/Gas.n,statsl[1].av()/Gas.n,T_av);


    if(write_distr)WriteDRR(distrfile);
    if(wr_film)WriteFilm(filmdir,Gas);





    if(fabs(1.-T_val/T_av)>cur_dc){
      if(cur_regime==SCALE){
         if(scale_vel && cur_iter<20){
           Gas.vel_scale(Gas.T/T_av);
           continue;
         }
      }
      if(cur_regime==TEST || cur_regime==SCALE || cur_regime==SWITCH_OFF){
         Gf1=Gf0;

         if(mc_equil){
           if(MC->momentum_depend)Gas.vel_scale(1.);
           MC->initial(Gas);
         }
         else Gas.ext_force=RandomForce;

         set_regime(ROUGH,dc1,rf_dt0,rf_nsteps);
         continue;
      }
    }// continues the cycle
    else{
      if(cur_regime==ACCEL){ // clearing meaningless statistics
          for(j=0;j<nstat;j++)stats[j]->clear();
          if(write_distr){
            DRRee.clear();
            DRRep.clear();
            if(non_symm)DRRpp.clear();
          }
      }

      if(cur_regime<=ROUGH){
         if(scale_vel)Gas.vel_scale(Gas.T/T_av);
         set_regime(SCALE,dc2,rf_dt0,tst_nsteps);
         continue;
      }

      if(cur_regime<=SCALE){
         if(e_stab){
           if(cur_iter<2 ||
	      fabs(Estb-statsl[2].av()/Gas.n)>fabs(Estb)*stab_acc){

             if(cur_iter>=2)printf("Total energy unstability %2.0f%%!\n",
				   100.*fabs(1.-statsl[2].av()/(Estb*Gas.n)));
             else printf("Energy stability check...\n");
	     Estb=statsl[2].av()/Gas.n;
             continue;
           }

         }
         if(rf_sw_off){
           Gf0=Gf1;
           set_regime(SWITCH_OFF,dc1,rf_dt0,sw_nsteps/sw_n);
           continue;
         }
         else if(e_stab){ // catching the average
           printf("Adjusting to average...\n");
           catch_average=1;
           Estab=statsl[2].av();
           MoveIt(t,wr_int,tst_nsteps,Gas,stats,statsl,nstat);
           if(catch_average){
             catch_average=0;
             printf("failed!\n");
             Estb=statsl[2].av()/Gas.n;
             continue;
           }
           else printf("OK\n");
         }

      }
      if(cur_regime==SWITCH_OFF){
			    if(cur_iter<sw_n)continue;
      }
      if(cur_regime<TEST && !no_test){
         Gas.ext_force=void_force1;
         set_regime(TEST,dc2,eq_dt,chk_nsteps);
         continue;
      }
      else break;
    }

   }while(iter<1000);

   if(iter>=1000){
     printf("Cannot reach equillibrium after 1000 cycles !\n");
     return -1.;
   }

   for(j=0;j<nstat;j++){
     *(stats[j])=statsl[j];
   }

   if(write_distr){
     DRRee.clear();
     DRRep.clear();
     if(non_symm)DRRpp.clear();
   }

   RFstatus=0.0;
   delete [] statsl;
   return t;
}



void WriteHeader(char *file, long wtype, Plasma &Gas,double Gamma, double T, char *dname,long dnstp){
 FILE *f1;
 f1=Err_fopen(file,"wb");
 char buff[400];
 sprintf(buff,"Trajectory data: %s\n, T: %f\nGamma: %f\nNparticles: %d\n",
     dname,T,Gamma,Gas.n);
 fwrite(buff,sizeof(char),200,f1);

 sprintf(buff,"%d %f %f %f %s %f %ld %f %ld\n",Gas.n,Gas.L,Gamma,T,dname,Gas.dt*dnstp,eq_nsteps/dnstp,Gas.Wpe,wtype);
 fwrite(buff,sizeof(char),200,f1);
 fwrite(Gas.x,sizeof(Vector_3),Gas.n,f1);
 fwrite(Gas.v,sizeof(Vector_3),Gas.n,f1);
 fclose(f1);
}




void WriteStep(char *file,long wtype, Plasma &Gas, long step){
 FILE *f1;
# ifndef UNIX
 f1=_fdopen(_open(file,O_BINARY|O_WRONLY),"wb");
# else
 f1=fdopen(open(file,O_WRONLY),"wb");
# endif
 if(!f1)serror("Can't reopen file %s",file);

 long dsz=gettypesize(wtype,Gas.n);

 fseek(f1,400+dsz*step,SEEK_SET);
 long pos=ftell(f1);
 if(pos!=400+dsz*step)printf("fseek failed!\n");
 ref_type=wtype|(VEL|COORD);

 new_trj=f1;
 anx=Gas.x;
 anv=Gas.v;
 npart=Gas.n;
 Refine_step();

 fclose(f1);

}

long CheckHeader(char *datafile,long wtype,Plasma &Gas,double Gamma,double T,char *ss, double *fdt, int stop, long &type){

  double l,g,t,im;
  int npart;

  char dname[100];

  npart=InitTrajectory(type,-1,datafile,dname,l,g,t,im);


  printf("Comparing:\n"
         "%15s%15s%15s\n","","Init","File");
  int cmp=0;

  printf("%15s%15s%15s\n","Dataname",dname,ss);
  //if(!strstr(dname,ss))cmp++;

  printf("%15s%15d%15d\n","Npart",Gas.n,npart);
  if(Gas.n !=npart)cmp++;

  printf("%15s%15f%15f\n","Gamma",Gamma,g);
  if(fabs(Gamma-g)>1e-5)cmp++;

  printf("%15s%15f%15f\n","T",T,t);
  if(fabs(T-t)>1e-5)cmp++;

  printf("%15s%15ld%15ld\n","Type",wtype,type);
  if(type!=0 && type!=wtype)cmp++;

  if(stop && cmp!=0)serror("Trajectory is inconsistent with Init data.\n");

  *fdt=dt;

  return nsteps;
}


void ReadGas(FILE *tFile,Plasma &Gas){
  fread(Gas.x,sizeof(Vector_3),Gas.n,tFile);
  fread(Gas.v,sizeof(Vector_3),Gas.n,tFile);
}

# define TEMPER  0x1
# define GAMMA   0x2
# define DENSITY 0x4

void AdjustGas(Plasma &Gas,Parameter *p, int np){
  int i;
  int  nad=0, adj_type=0;
  for(i=0;i<np;i++){
    if(strstr(p[i].s,"Gamma")){
      Gas.par_Gamma=get_pc(p[i]);
      adj_type|=GAMMA;
      nad++;
    }
    else if(strstr(p[i].s,"T")){
      Gas.par_T=get_pc(p[i]);
      adj_type|=TEMPER;
      nad++;
    }
    else if(strstr(p[i].s,"n")){
      Gas.par_density=get_pc(p[i]);
      adj_type|=DENSITY;
      nad++;
    }
  }
  if(nad<2)serror("Underdefined: not enough variable plasma parameters\n");
  if(nad>2)serror("Overdefined: too many contradictory variable parameters\n");

  switch(adj_type){
  case TEMPER|DENSITY:
    Gas.adjustTn(Gas.par_T,Gas.par_density);
    break;
  case TEMPER|GAMMA:
    Gas.adjustTG(Gas.par_T,Gas.par_Gamma);
    break;

  case GAMMA|DENSITY:
    Gas.adjustGn(Gas.par_Gamma,Gas.par_density);
    break;
  default:
    serror("Can not adjust plasma parameters.\n");
  }
}


void LoadFriedemann(char *ctrfile,Plasma &Gas){

  FILE *f1=Err_fopen(ctrfile,"rt");
  char str[800]; //,str1[200];
  int brk=0;
  do{

    if(feof(f1))serror("Can't find Friedemann's data in file %s.\n",ctrfile);
    fscanf(f1,"%s",str);
    if(str[0]=='#')fgets(str,800,f1);
    else brk=1;
  }while(!brk);

  int i,j;
  double val;
  for(i=0;i<12;i++)
    if(!fscanf(f1,"%lf",&val))serror("Bad Friedemann data\n");

  double me=0.9109534e-30;
  double unit_t=(1.602e-19)*(1.602e-19)/(4*M_PI*8.854e-12)*sqrt(me)*pow(1.38e-23*Gas.par_T*1e4,-3./2.);
  printf("unit_t: %e\n",unit_t);
  double ee;
  Vector_3 vav(0,0,0);

  int stat=0;
  do{
    ee=0.;
    vav=Vector_3(0,0,0);

    for(i=0;i<Gas.n;i++){
      for(j=0;j<3;j++){
       fscanf(f1,"%lf",&(Gas.x[i][j]));
       //printf("X %d, %d: %f",i,j,Gas.x[i][j]);
       Gas.x[i][j]*=5.2917655e-11/unit_l;
       Gas.x[i][j]-=Gas.L/2;
       //printf(" - %f %f\n",Gas.x[i][j],Gas.x[i][j]/Gas.L);
      }
      for(j=0;j<3;j++){
       fscanf(f1,"%lf",&(Gas.v[i][j]));
       //ee+=Gas.v[i][j]*Gas.v[i][j];
       //printf("V %d, %d: %f",i,j,Gas.v[i][j]);
       Gas.v[i][j]*=1.9939022e-24/unit_l/me*unit_t;
       ee+=Gas.v[i][j]*Gas.v[i][j];
       vav[j]+=Gas.v[i][j];
       //printf(" - %f\n",Gas.v[i][j]);
      }
      if(feof(f1)){
        serror("Bad Friedemann data: unexpected end-of-file\n");
        break;
      }
    }
    stat=1;
    for(i=0;i<13;i++){
      if(!fscanf(f1,"%lf",&val))break;
    }
  }while(i==13 && !feof(f1));

  vav/=Gas.n;
  printf("Read Friedem: Energy: %f, Vav: (%f, %f, %f) - %f\n", ee/3/Gas.n,vav[0],vav[1],vav[2],vav*vav/3);
  //for(i=0;i<Gas.n;i++)Gas.v[i]*=sqrt(3*Gas.n/ee);

  fclose(f1);
}

/*
int _matherr (struct exception *a)
{
  if (a->type == DOMAIN)
    if (!strcmp(a->name,"sqrt")) {
      a->retval = sqrt (-(a->arg1));
    return 1;
    }
  return 0;
} */




int main(int argc, char *argv[]){



# ifdef UNIX
 static char state[256];
 initstate(1997,state,256);
# else
  Exit_wait(1);
# endif

 Set_comment_char(';');

 StartTime();
 char *cfgfile="Params.dta";
 if(argc>1)cfgfile=argv[1];

 int no_remove=1;
 if(argc>2)if(strstr(argv[2],"new"))no_remove=0;

 Open_param_file(cfgfile);
 int nions,i;


 Read_param("Number of ions: %d",&nions);
 Read_param("Ionisation degree: %d",&i);
 ion_charge=i;

 double imass=1.;
 Set_stop(0);
 if(Read_param("Ion mass: %lf",&imass)){
   printf("Non-symmetric plasma specified!\n");
   non_symm=1;
 }

 char tmpstr[256];

 int sep_cm=0;
 if(Read_param("Center-of-mass for components: %s",tmpstr)){
   if(strstr(tmpstr,"separate")){
     printf("Plasma with separated center of masses specified!\n");
     sep_cm=1;
   }
 }



 Plasma *TheGas;  // allocation of the Gas
 double bdens, bTv;
 int bq;
 int bunch=Read_param("Bunch propagation: %lf, %lf, %d",&bdens,&bTv,&bq);
 if(bunch>0){
   printf("Bunch in plasma specified!\n");
   if(bunch!=3){
     msg_error("Invalid bunch specification!\n");
     exit(1);
   }
   bunch=1;
   bTv=sqrt(3*nions*i*bTv); // converting temperature to velocity
   TheGas=(Plasma *)new PlasmaBunch(bdens,bTv,bq,nions,i,imass);
 }
 else{
   TheGas=new Plasma(nions,i,imass);
   bunch=0;
 }
 Plasma &Gas=*TheGas;

 Gas.non_symm=non_symm;
 Gas.one_center=1-sep_cm;

 char dataname[50];
 Read_param("Data name: %s",dataname);
 strncpy(Gas.dataname,dataname,50);

 char ofile[256], pfile[256]="poten.dat";


 strcat(strcpy(pfile,dataname),".pot");

 Read_param("Output file: %s",ofile);
 if(strstr(ofile,"default")){
   strcpy(ofile,dataname);
   strcat(ofile,".eq");
 }

 Read_param("Log file: %s",logfile);
 if(strstr(logfile,"default")){
   strcpy(logfile,dataname);
   strcat(logfile,".log");
 }



 Parameter *p;
 int np=InitParameters(&p);


 double T,Gamma;

 potspec_t reader;
 reader.read_spec(cfgfile);
 Gas.potential=reader.potential;
 strncpy(Gas.charpot,reader.charpot,50);


 /*
 Read_param("Potential: %s",tmpstr);
 strncpy(Gas.charpot,tmpstr,50);
 if(strstr(tmpstr,"Kelbg"))Gas.potential=PotentialKELBG;
 else if(strstr(tmpstr,"Lennard-Johnes"))Gas.potential=PotentialJONES;
 else if(strstr(tmpstr,"Deutsch"))Gas.potential=PotentialDEUTSCH;
 else if(strstr(tmpstr,"Erf"))Gas.potential=PotentialERF;
 else if(strstr(tmpstr,"Cutoff")){
  if(!strcmp(tmpstr,"Cutoff1")){
    Gas.potential=PotentialCUT1;
  }
  else{
    Gas.potential=PotentialCUT;
    Read_param("Cutoff value*: %lf",&E_cut);
  }
 }
 else if(strstr(tmpstr,"ln")){
  Gas.potential=PotentialLN;
  Read_param("Cutoff value*: %lf",&E_cut);
 }
 else if(strstr(tmpstr,"table")){
  Gas.potential=PotentialTAB;
  Read_param("Potential table file: %s",tmpstr);
  Close_param_file();
  ReadPotential(tmpstr);
  Open_param_file(cfgfile);
 }
 else serror("Unknown potential type specified!\n");


 Read_param("Pauli part: %s",tmpstr);
 if(strstr(tmpstr,"n"))Pauli_part=0;


 double Lambda, Lambda_set;
 Read_param("R0: %s",tmpstr);
 if(strstr(tmpstr,"default"))Lambda_set=0.;
 else Lambda_set=atof(tmpstr);

 double Clam_ep=1.,Clam_ee=1;
 Read_param("e-e R0 coefficient: %lf",&Clam_ee);
 Read_param("e-p R0 coefficient: %lf",&Clam_ep);
 */

 int ask=0;
 Read_param("Dialog: %s",tmpstr);
 if(strstr(tmpstr,"y"))ask=1;

 auto_rf=0;
 Gf=10.;
 Read_param("Random force strength: %s",tmpstr);
 if(strstr(tmpstr,"auto"))auto_rf=1;
 else if(strstr(tmpstr,"fluct"))auto_rf=2;
 if(!sscanf(tmpstr,"%lf",&Gf) && auto_rf==0)serror("Can't read Random force strength\n");



 Read_param("Scale velocities: %s",tmpstr);
 if(strstr(tmpstr,"y"))scale_vel=1;

 Read_param("Delta: %lf",&delta);

 Read_param("Trajectory write interval: %ld",&wr_int);
 if(wr_int<=0)wr_int=-1;

 int new_rec=0;
 long wr_ions=0, wr_enseq=-1;
 Set_stop(0);
 if(Read_param("Ions write interval: %ld",&wr_ions)){
   new_rec=1;
   if(!Read_param("Electrons write sequence: %ld",&wr_enseq))wr_enseq=-1;
 }
 Set_stop(1);


 Read_param("Steps to check equillibrium: %ld",&chk_nsteps);
 Read_param("Steps with random force: %ld",&rf_nsteps);
 Read_param("Check steps with random force: %ld",&tst_nsteps);
 Read_param("Steps in equillibrium: %ld",&eq_nsteps);
 Read_param("Time step in equillibrium: %lf",&eq_dt);

 Read_param("Time step for random force: %lf",&rf_dt0);


 char trfile[256]="trajectory";

 int wr_tr=0;
 Set_stop(0);

 /*int pot_corr=0;
 if(Read_param("Potential correction: %s",tmpstr)){
   if(strstr(tmpstr,"y")){
     pot_corr=1;
     strcat(Gas.charpot," corr.");
   }                                       
 }*/


 if(Read_param("Total energy stability: %lf",&stab_acc))e_stab=1;


 if(Read_param("Positive cutoff: %lf",&E_negcut))neg_cut=1;
 else neg_cut=0;

 if(!Read_param("Random generator *:>",tmpstr))strcpy(tmpstr,"3");
 cList rndlist(tmpstr);

 int nrepeats=1;
 if(!Read_param("Repeats: %d",&nrepeats))nrepeats=1;


 char mdistrfile[256]="r-r.distrib";
 double rr_r0=0., rr_r1=-1.;
 if(Read_param("Write r-r distribution: %s",tmpstr)){
   if(strstr(tmpstr,"y")){
     write_distr=1;
     if(!Read_param("r-r file: %s",mdistrfile)||strstr(mdistrfile,"default")){
       strcpy(mdistrfile,"%s%d.rr");
     }
     if(Read_param("r-r range: %lf, %lf",&rr_r0,&rr_r1)!=2){
       rr_r0=0.;
       rr_r1=-1.;
     }
   }
 }


 if(Read_param("Soft step: %s",tmpstr)){
  if(strstr(tmpstr,"n"))soft_step=0;
 }

 if(Read_param("Soft random force: %s",tmpstr)){
  if(strstr(tmpstr,"y")){
    rf_sw_off=1;
    Set_stop(1);
    Read_param("Switch off steps: %ld",&sw_nsteps);
    Set_stop(0);
  }
  else rf_sw_off=0;
 }

 if(Read_param("Relative step: %s",tmpstr)){
  if(strstr(tmpstr,"y"))rel_step=1;
 }


 if(Read_param("Animation : %s",&tmpstr)){
   if(strstr(tmpstr,"y")){
      if(!Read_param("Film directory: %s",filmdir)||strstr(filmdir,"default")){
       strcpy(filmdir,"film/");
      }
   }
 }

 in_cs=-1.;
 Read_param("Initial cluster size: %lf",&in_cs);

 int restart=0,load_fried=0;
 int new_input=0;
 char inptrj[256];

 if(Read_param("Restart: %s",tmpstr)){
  if(strstr(tmpstr,"y")){
    restart=1;
    if(Read_param("Load Friedemann: %s",tmpstr)){
      if(strstr(tmpstr,"y"))load_fried=1;
    }
    if(Read_param("Input from: %s",inptrj))new_input=1;
  }

 }

 long wtype=0;
 if(Read_param("Trajectory file: %s",&trfile)){
  if(strstr(trfile,"default"))strcpy(trfile,"%s%d.r");



  Set_stop(1);
  wr_tr=1;

  Read_param("In output:>",tmpstr);
  if(strstr(tmpstr,"vel"))wtype|=VEL;
  if(strstr(tmpstr,"coord"))wtype|=COORD;
  if(strstr(tmpstr,"flow"))wtype|=FLOW;
  Set_stop(0);

 }
 else printf("Warning: no trajectory file\n");

 int mc_one=0;
 int mc_diff=0;
 int auto_adjust=1;
 double mc_inistep;

 int no_equilibr=0;
 Set_stop(1);
 if(Read_param("Equilibration procedure: %s",tmpstr)){
   if(strstr(tmpstr,"monte-carlo")){
     Set_stop(1);
     mc_equil=1;
     Read_param("One particle MC-step: %s",tmpstr);
     if(strstr(tmpstr,"y")){
       mc_one=1;
       // rf_dt0/=Gas.n;
     }
     Read_paramn(2,"MC step mode and value: %s %lf",tmpstr, &mc_inistep);
     if(strstr(tmpstr,"auto"))auto_adjust=1;
     else if(strstr(tmpstr,"stable"))auto_adjust=0;
     else serror("Unknown MC step mode.\n");
     Set_stop(0);
     mc_diff=0;
     if(Read_param("MC different temperatures: %s",tmpstr)){
       if(strstr(tmpstr,"y")){
	 if(mc_one)serror("Can use different MC temperatures\n"
			  "only in MC one-particle mode!\n");
	 mc_diff=1;
       }
     }
   }
   else{
     mc_equil=0;
     if(strstr(tmpstr,"off"))no_equilibr=1;
     else if(!strstr(tmpstr,"random-force")){
       serror("Unknown equilibration procedure: %s\n",tmpstr);
     }
   }
 }

 

 Set_stop(0);
 mc_calc=0;
 if(Read_param("Equilibrium calculation: %s",tmpstr)){
   if(strstr(tmpstr,"monte-carlo")){
     Set_stop(1);
     mc_calc=1;
     Read_param("One particle MC-step: %s",tmpstr);
     if(strstr(tmpstr,"y")){
       mc_one=1;
       // rf_dt0/=Gas.n;
     }
     Read_paramn(2,"MC step mode and value: %s %lf",tmpstr, &mc_inistep);
     if(strstr(tmpstr,"auto"))auto_adjust=1;
     else if(strstr(tmpstr,"stable"))auto_adjust=0;
     else serror("Unknown MC step mode.\n");
     Set_stop(0);
     mc_diff=0;
     if(Read_param("MC different temperatures: %s",tmpstr)){
       if(strstr(tmpstr,"y")){
	 if(mc_one)serror("Can use different MC temperatures\n"
			  "only in MC one-particle mode!\n");
	 mc_diff=1;
       }
     }
   }
 }

//# ifdef UNIX

 char out_dirs[250]="./",out_dirl[250]="./";
 if(Read_param("Data output directory: %s",out_dirs)){
   if(out_dirs[strlen(out_dirs)-1]!='/')strcat(out_dirs,"/");
   sprintf(tmpstr,out_dirs,dataname);
   strcpy(out_dirs,tmpstr);
# ifdef UNIX
   if(mkdir(out_dirs,S_IRWXU|S_IRGRP|S_IROTH)){
     if(errno!=EEXIST)serror("Can not create directory: %s.\n%s\n",
			     out_dirs,strerror(errno));
   }
   sprintf(tmpstr,"cp %s %s%s.cfg",cfgfile,out_dirs,dataname);
   //strcat(strcat(tmpstr,dataname),".cfg");
   if(system(tmpstr)==-1)
     printf("\nExec: %s\n",strerror(errno));
# else
   if(_mkdir(out_dirs)){
     if(errno!=EEXIST)serror("Can not create directory: %s.\n%s\n",
			     out_dirs,strerror(errno));
   }
   char cfgfilef[_MAX_PATH], out_dirsf[_MAX_PATH];
   _fullpath(cfgfilef, cfgfile, _MAX_PATH);
   _fullpath(out_dirsf, out_dirs, _MAX_PATH);
   sprintf(tmpstr,"copy %s %s%s.cfg",cfgfilef,out_dirsf,dataname);
   //strcat(strcat(tmpstr,dataname),".cfg");
   if(system(tmpstr)==-1)
     printf("\nExec: %s\n",strerror(errno));
# endif

   


   strcpy(ofile,strcat(strcpy(tmpstr,out_dirs),ofile));
   strcpy(mdistrfile,strcat(strcpy(tmpstr,out_dirs),mdistrfile));
   strcpy(sfile,strcat(strcpy(tmpstr,out_dirs),sfile));
   strcpy(pfile,strcat(strcpy(tmpstr,out_dirs),pfile));
   if(wr_film){
     strcpy(filmdir,strcat(strcpy(tmpstr,out_dirs),ofile));
# ifdef UNIX
     if(mkdir(filmdir,S_IRWXU|S_IRGRP|S_IROTH)){
# else
     if(_mkdir(filmdir)){
# endif
       if(errno!=EEXIST)serror("Can not create directory: %s.\n%s\n",
			       filmdir,strerror(errno) );
     }
     strcat(filmdir,dataname);
   }
 }
 if(Read_param("Process output directory: %s",out_dirl)){
   if(out_dirl[strlen(out_dirl)]!='/')strcat(out_dirl,"/");
   sprintf(tmpstr,out_dirl,dataname);
   strcpy(out_dirl,tmpstr);


# ifdef UNIX
   if(mkdir(out_dirl,S_IRWXU|S_IRGRP|S_IROTH)){
# else
   if(_mkdir(out_dirl)){
# endif 
     if(errno!=EEXIST)serror("Can not create directory: %s.\n%s\n",
			     out_dirl,strerror(errno));
   }


   strcpy(logfile,strcat(strcpy(tmpstr,out_dirl),logfile));
   strcpy(trfile,strcat(strcpy(tmpstr,out_dirl),trfile));
 }

//# endif  // UNIX

 int spwn_trj=0;

 if(Read_param("Spawn trajectories: %s",tmpstr)){
   if(strstr(tmpstr,"y")){
    spwn_trj=1;
    //Read_paramn(1,"Spawn frame: %d",&spwn_frame);
    no_test=1; // no RF test for equilibrartion
   }
 }

 Gas.stable_ions=0;
 if(Read_param("Stable ions: %s",tmpstr)){
   if(strstr(tmpstr,"y")){
    Gas.stable_ions=1;

   }
 }

 int repc_i=0;
 if(!Read_param("Repeat counter start: %d",&repc_i))repc_i=0;

 int irescale=0;
 if(Read_param("Ion velocity rescale: %s",tmpstr)){
  if(strstr(tmpstr,"y"))irescale=1;
  if(strstr(tmpstr,"reset"))irescale=2;
 }

 /*
 int limrescale=0;
 if(Read_param("Rescale on electron temperature reached: %lf  %ld",&limTe,&limstpe)==2){
  limrescale=1;
  limspec=0x2;
  limrescale_switch(0);
 } */

 int limrescale=0;
 int inc_mes=0;

 if(Read_param("Incremental measurement (T0,dT,mes_steps): %lf,%lf,%ld",&incT0,&incdT,&incStp)==3){
   inc_mes=1;
   fixT=incT0;
   limstpe=1;
   limstpi=1;
 }



 Set_stop(1);



 Gas.ini_Te=Gas.ini_Ti=1.;  // by default equal temperatures
 Read_param("Initial velocity distribution: %s",tmpstr);
 if(strstr(tmpstr,"maxwell"))in_distr=MAXWELL;
 else if(strstr(tmpstr,"max_polak"))in_distr=MAXWELL_P;
 else if(strstr(tmpstr,"zero")){
   if(mc_equil){
     in_distr=MAXWELL;
     printf("Warning: setting 'maxwell' initial vel. distribution for MC!\n");
   }
   else in_distr=ZEROVEL;
 }
 else if(strstr(tmpstr,"separate")){
   in_distr=SEPARATE;
   int &ndistr=Gas.idistr;
   int nr;
   char relstr[200];
   for(i=0;i<2;i++){
    if(i==0)nr=  Read_param("Electron velocity distribution: %s  %lf  %s",tmpstr,&Gas.ini_Te,relstr);
    else{
     nr= Read_param("Ion velocity distribution: %s  %lf %s",tmpstr,&Gas.ini_Ti,relstr);
     Gas.edistr=ndistr;
    }
    if(nr<2){
      serror("Can't read velocity distribution parameters!\n");
    }
    if(nr>=2){
      if(strstr(relstr,"abs")){
        if(i==0)Gas.rel_Te=0;
        else Gas.rel_Ti=0;
      }
    }

    if(strstr(tmpstr,"maxwell"))ndistr=MAXWELL;
    else if(strstr(tmpstr,"max_polak"))ndistr=MAXWELL_P;
    else if(strstr(tmpstr,"zero"))ndistr=ZEROVEL;
    else {
     printf("Warning: unknown distribution '%s', setting to 'zero'\n",tmpstr);
     ndistr=ZEROVEL;
    }
   }
 }
 else{
   printf("Warning: unknown distribution '%s', setting to 'zero'\n",tmpstr);
   in_distr=ZEROVEL;
 }

 

 Close_param_file();


 Statistics sT(&Gas.T),sEcoul(&Gas.Ecoul),sEpotent(&Gas.Epotent),sQuant(&Gas.Quant), sEtot(&Gas.Etot);
 Statistics sTi(&Gas.Ti),sTe(&Gas.Te);
 const int nstat=7;
 Statistics *stats[nstat]={&sT,&sEcoul,&sEpotent,&sQuant,&sEtot,&sTi,&sTe};

 SetTableform(GNU);
 if(no_remove){
   no_remove=CheckData(sfile,ofile,dataname,np,p,ask);
 }

 if(restart && !wr_tr && !new_input)serror("No trajectory file specified, cannot restart.\n");
 if(restart && new_input && wr_tr)
   if(!strcmp(trfile,inptrj))
     serror("Equal names for input and output trj-files.\n");


 FILE *f1;
 if(!no_remove || no_remove==2){
  if(!no_remove)f1=Err_fopen(ofile,"wt");
  else f1=Err_fopen(ofile,"at");
  BeginFrame(f1,dataname,np,p);
  fprintf(f1,"%9.9s %9.9s %9.9s %9.9s %9.9s %9.9s %9.9s %9.9s\n",
             "T","Ecoul","Epotent","Equant","dT","dEcoul","dEpotent","dEquant");
  fclose(f1);
 }

 double t=0;

 f1=fopen(logfile,"rt");
 if(f1){
   if(!restart && !no_remove){
     fclose(f1);
     remove(logfile);
   }
   else {// reading 'log' time
     fseek(f1,0,SEEK_END);
     long pos=ftell(f1)-4;
     do{
       fseek(f1,pos,SEEK_SET);
       if(fgetc(f1)=='\n'){
	        fscanf(f1,"%lf",&t);
         //printf("new t: %f\n",t);
	        break;
       }
       pos--;
     }while(pos>=0);
     fclose(f1);
   }
 }


 int ccount=CurrentCount(np,p);
 WriteStatus(!ccount,sfile,dataname,np,p);



 do{

   //Gamma=get_pc(p[0]);
   //T=get_pc(p[1]);



  char ctrfile[256];
  sprintf(ctrfile,trfile,dataname,ccount);
  char ss[100];
  sprintf(ss,"%s%d",dataname,ccount);
  sprintf(distrfile,mdistrfile,dataname,ccount);

  // StatusLine(str,np,p);
  //show_status(440,str);

  long i,n,m,nf=0;

  //Gas.adjustTG(T,Gamma);

  AdjustGas(Gas,p,np);
  T=Gas.par_T;
  Gamma=Gas.par_Gamma;

  double me=0.9109534e-30;
  double qe=1.602e-19;
  double Kb=1.38e-23;
  double unit_t=qe*qe/(4*M_PI*8.854e-12)*sqrt(me);


  // overcoming g++ bug with -O3
  double jojo=Kb*T*1.e4;
  //printf("c2: %e\n",jojo);
  jojo=jojo*jojo*jojo;
  jojo=sqrt(jojo);
  //printf("c21: %e\n",jojo);
  unit_t/=jojo;
  //printf("c3:\n");

  double unit_h=Kb*T*1e4*unit_t;
  double h_qwer=1.0545887e-34;

  reader.calc_lambda(T,Gas.ini_Te,Gas.ini_Ti,Gas.mass);

  /*
  if(!non_symm){
   if(Lambda_set!=0.0){
     Lambda=Lambda_set/unit_l;
     equal_lambda=0;
   }
   else {
     Lambda=0.179*sqrt(T);
     equal_lambda=1;
   }
   Lambda_pauli=0.179*sqrt(T);
   //non_symm=0;
   Lambda_ee=Lambda_pp=Lambda*Clam_ee;
   if(fabs(Clam_ee-1.)>1e-5)equal_lambda=0;
   Lambda_ep=Lambda*Clam_ep;
  }
  else{
   if(Lambda_set!=0.0){
     serror("Can not setup lambda for nonsymmetric plasma.\n");
     //Lambda=Lambda_set/unit_l;
     //Lambda_ep=Lambda*Clam_ep;
     //Lambda_ee=Lambda*Clam_ee;
     //equal_lambda=0;
   }
   else{
     if(strstr(Gas.charpot,"Deutsch")){
       printf("Setting lambda values for Deutsch potential !!!\n");
       Lambda_pauli=h_qwer/(unit_h);
       Lambda=Lambda_pauli/sqrt(2*M_PI);
       double m_ee=0.5, m_pp=0.5*Gas.mass, m_ep=Gas.mass/(1.+Gas.mass);

       Lambda_ee=Lambda/sqrt(m_ee);
       Lambda_pp=Lambda/sqrt(m_pp);
       Lambda_ep=Lambda/sqrt(m_ep);
     }
     else { //if(strstr(Gas.charpot,"Kelbg")){
       printf("Setting lambda values for Kelbg potential !!!\n");

       Lambda_pauli=h_qwer/(unit_h);
       Lambda=Lambda_pauli;
       double m_ee=0.5, m_pp=0.5*Gas.mass, m_ep=Gas.mass/(1.+Gas.mass);
       double t_ep=Gas.ini_Te;
       double t_pp=Gas.ini_Ti;
       double t_ee=Gas.ini_Te;

       Lambda_ee=Lambda/sqrt(2*m_ee*t_ee);
       Lambda_pp=Lambda/sqrt(2*m_pp*t_pp);
       Lambda_ep=Lambda/sqrt(2*m_ep*t_ep);
     }

   }
  } */


  if(rel_step)rf_dtcoeff=1./Gas.Wpe;

  double kk=(1.602e-19)*(1.602e-19)/(4*M_PI*1.38e-23*T*1e4*8.854e-12);



  printf("Simulation parameters:\n");
  printf("Gamma =%f\n"
	 "L=%f=%12e m\n"
	 "lambda=%f=%12e m\n"
	 "1/Wp=%f\n"
	 "Rd=%f\n",Gamma,Gas.L, Gas.L*kk,reader.Lambda,reader.Lambda*kk,1./Gas.Wpe,RDebye);
  printf("Density = %1.2e cm^(-3)\n"
	 "T= %f K\n",Gas.par_density*1e19,T*1e4);

  printf("Time step relations:\n");
  printf("dtrf*Wp= %f, dteq*Wp= %f\n", (rel_step ? rf_dt0 : rf_dt0*Gas.Wpe),
	 (rel_step ? eq_dt : eq_dt*Gas.Wpe));


  double t_inter=RDebye/sqrt(2*getEmax(Gas));


  printf("dtrf/tint= %f, dteq/tint= %f\n", rf_dt0*rf_dtcoeff/t_inter,
	 eq_dt*rf_dtcoeff/t_inter);


  reader.calc_correction(Gas.par_T);
  reader.write_pot(pfile, Gas.L);

  /*
  //if(ccount==1)t/=Gas.Wpe;
  if(pot_corr){
    //Correction(IONION,Gas.potential,Gas.par_T);
    Correction(IONELC,Gas.potential,Gas.par_T);
    Correction(ELCELC,Gas.potential,Gas.par_T);
  }*/

  Statistics rs[nstat];



  int repi=0;
  int repc=repc_i;

  do{ //through nrepeats

    printf("\nStarting calculation #%d...\n",repc);
    flm_count=0;

    if(repc>0){
      repi++;
      sprintf(tmpstr,"%02d",repc);
      if(repi>1){
       distrfile[strlen(distrfile)-2]=0;
       ctrfile[strlen(ctrfile)-2]=0;
      }
      strcat(distrfile,tmpstr);
      strcat(ctrfile,tmpstr);
    }


    if(write_distr){
      DRRee.init(rr_r0,(rr_r1>0 ? rr_r1 : Gas.L/2),400);
      DRRep.init(rr_r0,(rr_r1>0 ? rr_r1 : Gas.L/2),400);
      if(non_symm)DRRpp.init(0,(rr_r1>0 ? rr_r1 : Gas.L/2),400);
    }

    long ftype=wtype;


    if(restart){
      printf("restarting...\n");
      if((mc_equil && spwn_trj) || mc_calc){
        mc_equil=-1;
        MC = new mc_simm(Gas,mc_inistep*Gas.L /*Gas.L/Gas.n*/,0.5,
			 mc_one,mc_diff);
        if(!MC)serror("Cannot allocate MCsimm\n");
        MC->auto_adjust=auto_adjust;
      }
      else mc_equil=0;

      if(!new_input)strcpy(inptrj,ctrfile);

      if(load_fried){
        LoadFriedemann(inptrj,Gas);
        Gas.dt=eq_dt;
        if(rel_step)Gas.dt/=Gas.Wpe;

        if(eq_nsteps<wr_int)eq_nsteps=wr_int;
	m=wr_int;
	n=eq_nsteps/m;
	Gas.ext_force=void_force1;
	wr_tr=0;
	ftype=0;
      }
      else{
	Gas.dt=(rel_step ? eq_dt/Gas.Wpe : eq_dt);
	Trajectory.Check(inptrj,wtype,Gas,wr_int,wr_ions,wr_enseq,new_input);
	if(Trajectory.wtype !=0 && !new_input){
	  Gas.dt=Trajectory.AdjustInterval(Gas.dt);
	}

	long stp;
	if(!new_input){
	  stp=Trajectory.ReloadGas(Gas,1);
	  t=stp*Trajectory.file_dt();
	  nf=(long)(t/Gas.dt/wr_int+0.01);
	}
	else{
	  stp=Trajectory.ReloadGas(Gas,0);
	  t=0;
	  nf=0;
	}
	printf("t= %f, nf= %ld,  %f\n",t,nf,(t/Gas.dt/wr_int));

	//serror("Restart is not yet implemented !\n");
        Gas.ext_force=void_force1;
      }
      //restart=0;
    }
    else{
      rand_init=rndlist.step();
      if(rand_init<0){
        rndlist.rewind();
        rand_init=rndlist.step();
      }

      strcat(distrfile,"e");


      if((mc_equil && (!spwn_trj || (repc==repc_i && spwn_trj))) || mc_calc){
        mc_equil=1;
        MC = new mc_simm(Gas,mc_inistep*Gas.L /*Gas.L/Gas.n*/,0.5,
			 mc_one,mc_diff);
        if(!MC)serror("Cannot allocate MCsimm\n");
        MC->auto_adjust=auto_adjust;
      }

      if(!no_equilibr){
        if(!spwn_trj || (spwn_trj && repc==repc_i)){
          limrescale_switch(0);
          Equillibrium(t,Gas,stats,nstat);
          limrescale_switch(limrescale);
        }
      }
      else{
        Gas.r0=0.05;
        Gas.init_config(in_cs,in_distr);
      }

      distrfile[strlen(distrfile)-1]=0; // deleting 'e'

      if(mc_equil){
       if(spwn_trj){
        Gas.init_vel(in_distr);
       }
       else if(!mc_calc){
	 delete MC;
	 mc_equil=-1;
       }
       else
         mc_equil=1;
      }

      Gas.dt=eq_dt;
      if(rel_step)Gas.dt/=Gas.Wpe;

    }
    if(irescale){
     if(irescale==1)Gas.ivel_scale(1.);
     else Gas.init_vel(in_distr);
    }

    if(eq_nsteps<wr_int)eq_nsteps=wr_int;
    if(wr_int>0){
      n=eq_nsteps/wr_int;
      m=wr_int;
    }
    else if(eq_nsteps>=500){
      n=eq_nsteps/500;
      m=500;
    }
    else{
      n=1;
      m=eq_nsteps;
    }

    if(wr_tr && (!restart || new_input || (restart && ftype==0))){

      //WriteHeader(ctrfile,wtype,Gas,Gamma,T, ss,m);

      // initializing PlasmaRec
      Trajectory.Clear();
      Trajectory.Init(ctrfile,wtype,Gas,wr_int,wr_ions,wr_enseq);
    }
    if(wr_tr && (restart && ftype==0)){
      //WriteStep(ctrfile,wtype,Gas,0);
    }

    Statistics statsl[nstat];

    Trajectory.valid=wr_tr;
    // new cycle

    if(bunch){
      PlasmaBunch *pb=(PlasmaBunch *)&Gas;
      pb->start_bunch();
    }
    for(i=nf;i<n;i++){
      double term_c=1.;
      if(Gas.stable_ions)term_c=(double)Gas.ne/Gas.n;
      printf("Equilibrium calculation: %f%% complete, T= %f Et=%f\n",(i+1)*100./n,Gas.T,Gas.T*3/2*term_c+Gas.Epotent/Gas.n);

      if(StopStatus(sfile,3)){
       	StopStatus(sfile,-1);
        serror("Program interrupted!\n");
      }

      if(mc_equil && spwn_trj){
        mc_equil=1;
        double ratio=MC->get_ratio();
        if(ratio<1e-10)printf("Estimated randomization steps: infinity\n");
        else printf("Estimated randomization steps: %d\n",(int)(Gas.n*Gas.L/MC->dc[1]/ratio));
      }
      MoveIt(t,m,m,Gas,stats,statsl,nstat);


      //if(wr_tr)WriteStep(ctrfile,wtype,Gas,i);
      if(write_distr)WriteDRR(distrfile);
      if(wr_film)WriteFilm(filmdir,Gas);
    }
    if(bunch){
      PlasmaBunch *pb=(PlasmaBunch *)&Gas;
      pb->stop_bunch();
    }

    if(write_distr)WriteDRR(distrfile);

    for(i=0;i<nstat;i++)rs[i]+=*(stats[i]);
    repc++;

    if(new_input)restart=0;
  }while(repc<nrepeats);

  if(spwn_trj)delete MC;

  f1=Err_fopen(ofile,"at");
  MiddleFrame(f1,np,p);
  int gn=Gas.n;
  //fprintf(f1,"%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
  //       sT.av()*T,sEcoul.av()/gn,sEpotent.av()/gn,sQuant.av()/gn,
  //       sT.dev()*T,sEcoul.dev()/gn,sEpotent.dev()/gn,sQuant.dev()/gn);

  fprintf(f1,"%9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",
         rs[0].av()*T,rs[1].av()/gn,rs[2].av()/gn,rs[3].av()/gn,
         rs[0].dev()*T,rs[1].dev()/gn,rs[2].dev()/gn,rs[3].dev()/gn);

  fclose(f1);

  ccount=CycleCount(np,p);
  WriteStatus(!ccount,sfile,dataname,np,p);

  restart=0;
  repc=0;
 }while(ccount);
 return 0;
}






























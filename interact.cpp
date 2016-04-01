#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "common.h"
#include "rparam.h"

#include "interact.h"

# ifndef UNIX
# include "erf.h"
# endif


// potentials must return k=-dE/dr


double Lambda_pp,Lambda_ee,Lambda_ep,Lambda_pauli;

int Pauli_part=1, equal_lambda=1, non_symm=0;

double ion_charge=1.;

int neg_cut=0;
int corrected[4]={0,0,0,0};

double corr_value[4]={0.,0.,0.,0.};
double corr_valueq[4]={0.,0.,0.,0.};

double E_cut=1.;
double E_negcut=1.;
double E_jones=1.;

double gkR0_ee=1, gkaR0_ee=1;
double gkR0_ei=1, gkaR0_ei=1;
double gkR0_ii=1, gkaR0_ii=1;
int gkpot_ee=GK_COUL, gkpot_ei=GK_GK, gkpot_ii=GK_KELBG;

double const SQPI=1.772453850905516, SQ2=1.414213562373095;



double integ(double xi){
  double dt=1e-3;
  double t=dt,sum=0.,ds,r;
  do{
   if(xi>0){
    ds=t*exp(-t*t)/(1-exp(-xi*M_PI/t))*dt;
   }
   else{
    r=exp(xi*M_PI/t);
    ds=t*exp(-t*t)*r/(r-1.)*dt;
   }

   sum+=ds;
   t+=dt;
  }while(t<20.);

  return 4*sqrt(M_PI)*xi*sum;
}


double sum(double xi){
  if(xi<0)return 0.;
  double s=0;
  int i;
  for(i=1;i<100;i++){
    s+=(1.+xi*xi/(4*i*i))/(i*i*i);
  }
  return (double)s*sqrt(M_PI)*xi*xi*xi;
}

double sum_adv(double xi){
  if(xi<0)return 0.;
  double s=0,ds;
  double i;
  for(i=1;i<1000.;i+=1){
    ds=(1.+xi*xi/(4*i*i))/(i*i*i);
    s+=ds;
    if(ds<1e-4)break;
  }
  if(i>=1000.)
    printf("! Sum accuracy is: %f\n",ds);

  return (double)s*sqrt(M_PI)*xi*xi*xi;
}




double fcorr(int type,double T){
  double f;
  double me=0.9109534e-30;
  double kb=1.389622e-23;
  double hq=1.0545887e-34;

  double t_k=(1.602e-19)*(1.602e-19)/(4*M_PI*8.854e-12)*sqrt(me)*pow(kb*1e4,-3./2.);
  double h_k=hq/(1e4*kb*t_k);
  //printf("h_k: %f\n",h_k);
  double lambda_k=h_k; ///sqrt(2.);
  double xi_k=1./lambda_k;
  
  double xi=xi_k/sqrt(T);


  if(type==IONELC || type==ELCION){
    f=sum_adv(xi)+integ(xi);
  }
  else {
    xi=-xi;
    f=0.5*integ(xi);
  }
  
  if(f<1e-32)f= -32;
  else if(f>1e32)f=32;
  else f=log(f);

  return f;
}


double Correction(int type,double potential(int,double,double&,double&),double T){

  double f=fcorr(type,T);
  double dT=1e-5;
  double fder;

  fder=(fcorr(type,T+dT)-fcorr(type,T-dT))/(2*dT);
  double dEpot, dEquant;
  
  corrected[type]=0;
  potential(type,1e-6,dEpot,dEquant);
 
  if(type==IONELC || type==ELCION){
    corr_value[IONELC]=corr_value[ELCION]= -f- dEpot;
    corr_valueq[IONELC]=-f+T*fder-dEquant;
    corrected[IONELC]=corrected[ELCION]=1;
  } 
  else{
    corr_value[type]= -f- dEpot;
    corr_valueq[type]=-f+T*fder-dEquant;
    corrected[type]=1;
    if(!non_symm){
      corr_value[3-type]= -f- dEpot;
      corr_valueq[type]=-f+T*fder-dEquant;
      corrected[3-type]=1;
    }
  }

 
  return f;
}




  



void WritePotential(const char *file,double start, double L,int n_res,double  f(int,double ,double&,double&)){
 FILE *f1=Err_fopen(file,"wt");
 double F;
 double Ep,Eq,x;

 fprintf(f1,"# Potential data\n");

 if(!non_symm)fprintf(f1,"# Symmetric: yes\n");
 else fprintf(f1,"# Symmetric: no\n");
 fprintf(f1,"#1-r 2-U(r) 3-force(r) 4-Uquant(r) 5-force_num(r)\n");

 int type=ELCELC;
 
 double dx=(L-start)/(n_res+1);
 double dxn=1e-4;
 if(dx/2<dxn)dxn=dx/2;
 for(int i=0;i<3;i++){
   if(i==0){
     type=ELCELC;
     fprintf(f1,"# e-e interaction:\n");
   }
   else if(i==1){
     if(!non_symm)
       continue;
     type=IONION;
     fprintf(f1,"\n\n# i-i interaction:\n");
   }
   else {
     type=ELCION;
     fprintf(f1,"\n\n# e-i interaction:\n");
   }

   for(x=start+dx;x<=L;x+=dx){
     double Fn;
     f(type,x+dxn,Fn,Eq);
     f(type,x-dxn,F,Eq);
     Fn=(Fn-F)/(2*dxn);
     F=f(type,x,Ep,Eq);
   
     fprintf(f1,"%e %e %e %e %e\n",x,Ep,F,Eq,-Fn);
   }
 }
 fclose(f1);
} 


double PotentialKELBG(int type,double R,double& dEpotent, double &dQuant){  
  double R_inv=1/R; 
  double R_inv2=R_inv*R_inv;
  double k;
  double Rlam,Rlam2,ex,t1;
   
  switch(type){

  case IONION:
  case ELCELC:
   

    if(type==ELCELC)Rlam=R/Lambda_ee;
    else Rlam=R/Lambda_pp;
 
    Rlam2=Rlam*Rlam;
    ex=exp(-Rlam2);
    t1=SQPI*Rlam*erfc(Rlam);

    dEpotent=R_inv*(1.-ex+t1);
    

    if(Pauli_part && !(non_symm && type==IONION)){ //only for electrons or protons (not for ions)
      if(!equal_lambda){
	double Rlam_p=R/Lambda_pauli;
	double Rlam2_p=Rlam*Rlam;	
	double ex_p=exp(-Rlam2);
      	
	dEpotent+=0.5*ex;
	k=R_inv2*(1.-ex)+ex_p*Rlam_p/Lambda_pauli;
	t1+=-Rlam2_p*ex_p;  
      } 
      else{
	dEpotent+=0.5*ex;
	k=R_inv2*(1.-(1.-Rlam2*R)*ex);
	t1+=-Rlam2*ex;
      }
    }
    else k=R_inv2*(1.-ex);
    dQuant=dEpotent-0.5*t1*R_inv;

    break;

  default:
  case IONELC:
  case ELCION:
    Rlam=R/Lambda_ep;
    Rlam2=Rlam*Rlam;
    ex=exp(-Rlam2);
    t1=SQPI*Rlam*erfc(Rlam);


    dEpotent=-R_inv*(1.-ex+t1);
    t1=-t1;
    k=-R_inv2*(1.-ex);
   
    dQuant=dEpotent-0.5*t1*R_inv;
    break;
  }
  //dQuant=R_inv-0.5*t1;
  if(corrected[type]){
    dEpotent+=corr_value[type]*ex;
    t1=2*Rlam2*corr_value[type]*ex;
    k+=t1*R_inv;
    dQuant+=ex*corr_valueq[type] - 0.5*t1;
  }
  
  return k;
}



double PotentialERF(int type,double R,double& dEpotent, double &dQuant){
  double R_inv=1/R;
  double R_inv2=R_inv*R_inv;
  double k=0;
  double Rlam,Rlam2,ex,t1,erff;
  
  
  switch(type){

  case IONION:
  case ELCELC:
    if(type==ELCELC)Rlam=R/Lambda_ee;
    else Rlam=R/Lambda_pp;

    Rlam2=Rlam*Rlam;
    ex=exp(-Rlam2);
    t1=SQ2*sqrt(ex)/(SQPI*Lambda_ee);
    erff=erf(Rlam/SQ2);

   
    dEpotent=R_inv*erff; 
    k=R_inv2*erff-R_inv*t1;
    if(Pauli_part && !(non_symm && type==IONION)){ //only for electrons or protons (not for ions) 
      if(!equal_lambda){
	Rlam=R/Lambda_pauli;
	Rlam2=Rlam*Rlam;	
	ex=exp(-Rlam2);
      }	
      dEpotent+=0.5*ex;
      k+=ex*Rlam/Lambda_pauli;
      t1+=-Rlam2*ex;  
    }    
    dQuant=dEpotent-0.5*t1;
    break;

  case IONELC:
  case ELCION:
    Rlam=R/Lambda_ep;
    Rlam2=Rlam*Rlam;
    ex=exp(-Rlam2);
    t1=SQ2*sqrt(ex)/(SQPI*Lambda_ep);
    erff=erf(Rlam/SQ2);

    
    dEpotent=-R_inv*erff;
    t1=-t1;
    k=-R_inv2*erff+R_inv*t1;
    dQuant=dEpotent-0.5*t1;
    break;
  }
  //dQuant=dEpotent-0.5*t1;
  return k;
}



double PotentialCUT(int type,double R,double& dEpotent, double &dQuant){
  double R_inv=1/R; 
  double k=0.;
  
  switch(type){

  case IONION:
    dEpotent=R_inv*ion_charge*ion_charge;
    k=dEpotent*R_inv;
    dQuant=0;
    break;
  case ELCELC:

    if(neg_cut && R_inv>E_negcut){
      dQuant=dEpotent=E_negcut;
      k=0.;
    }
    else{
      dEpotent=R_inv;
      k=R_inv*R_inv;
      dQuant=0;
    }
    if(Pauli_part && !(non_symm && type==IONION)){ //only for electrons or protons (not for ions)
      double Rlam=R/Lambda_pauli;
      double ex=exp(-Rlam*Rlam);
      dEpotent+=0.5*ex;
      k+=ex*Rlam/Lambda_pauli;
      dQuant+=-0.5*Rlam*Rlam*ex;
    }

    break;
  case IONELC:
  case ELCION:
    if(R_inv>E_cut){
      dQuant=dEpotent=-E_cut*ion_charge;
      k=0.;
    }
    else{
      dEpotent=-R_inv*ion_charge;
      k=dEpotent*R_inv;
      dQuant=0;
    }
    break;
  }
  dQuant=dEpotent-dQuant;
  return k;
}


double PotentialLN(int type,double R,double& dEpotent, double &dQuant){
  double R_inv=1/R;
  double k=0.;

  switch(type){

  case IONION:
  case ELCELC:
    if(neg_cut && dEpotent > E_negcut){
      dQuant=dEpotent=E_negcut;
      k=0.;
    }
    else {
      dEpotent=R_inv;
      k=R_inv*R_inv;
      dQuant=0;
    }

    if(Pauli_part && !(non_symm && type==IONION)){ //only for electrons or protons (not for ions)
      double Rlam=R/Lambda_pauli;
      double ex=exp(-Rlam*Rlam);
      dEpotent+=0.5*ex;
      k+=ex*Rlam/Lambda_pauli;
      dQuant+=-0.5*Rlam*Rlam*ex;
    }
    dQuant=R_inv-dQuant;
    break;
  case IONELC:
  case ELCION:

    int no_force=0;
    if(R_inv > E_cut){
      R_inv=E_cut;
      no_force=1;

    }

    double ex=exp(-R_inv);
    double sqR_inv=sqrt(R_inv);
    double erfcf=erfc(sqR_inv);
    double R32=R_inv*sqR_inv;
    double R52=R32*R_inv;
    double corr=(2.*sqR_inv+4./3.*R32+8./15.*R52)/SQPI;
    double inn=erfcf+ex*corr;

    dEpotent=-R_inv-log(inn);
    double t1=(ex/inn)*(-8./15.*R52*R_inv*R_inv)/SQPI;

    if(no_force){
      dQuant=dEpotent;
      k=0;
    }
    else{
      k=-R_inv*R_inv-t1;
      dQuant=t1*R;
    }
    dQuant=-R_inv-dQuant;
    break;
  }
  //dQuant=-R_inv-dQuant;
  return k;
}

double PotentialDEUTSCH(int type,double R,double& dEpotent, double &dQuant){
  double R_inv=1./R;
  double R_invsq=R_inv*R_inv;

  double k=0.;

  switch(type){

  case IONION:
   {
    double Rlam=R/Lambda_pp;
    double e=exp(-Rlam);
    dEpotent=R_inv*(1.-e);
    k=R_invsq-(R_invsq+R_inv/Lambda_pp)*e;
   }
  break;
  case ELCELC:
   {
    double Rlam=R/Lambda_ee;
    double e=exp(-Rlam);
    dEpotent=R_inv*(1.-e);
    k=R_invsq-(R_invsq+R_inv/Lambda_ee)*e;
   }
  break;

  case IONELC:
  case ELCION:
   {
    double Rlam=R/Lambda_ep;
    double e=exp(-Rlam);
    dEpotent=-R_inv*(1.-e);
    k=-R_invsq+(R_invsq+R_inv/Lambda_ep)*e;
   }
   break;
  }

  dQuant=0.;
  return k;
}


inline double PotGK(int type, double Zz, int gkpot, double gkR0, double gkaR0, double R,double& dEpotent, double &dQuant){
  if(gkpot==GK_GK){
    // e-e is according to Gurski, Krasko FTT, 11, 3016 (1969) 
    double R_inv=1./R;
    double e=exp(-R/gkR0);
    double er=e*R_inv;
    double de1=R_inv-er;
    double de2=gkaR0*e;
    double aZz=fabs(Zz);
    dEpotent=de1*Zz+de2*aZz;
    dQuant=0.;
    return de1*Zz*R_inv+(-er*Zz+de2*aZz)/gkR0;
  }
  else if(gkpot==GK_KELBG){
    return PotentialKELBG(type,R,dEpotent, dQuant);
  }
  // else COULOMB
  double R_inv=1./R;
  double R_invsq=R_inv*R_inv;
  dEpotent=R_inv*Zz;
  dQuant=0.;
  return R_invsq*Zz;
}


double PotentialGURKR(int type,double R,double& dEpotent, double &dQuant){
  if(type==ELCELC) // e-e interaction 
    return PotGK(type,1.,gkpot_ee,gkR0_ee,gkaR0_ee,R,dEpotent,dQuant);
  if(type==IONION)  // i-i interaction 
    return PotGK(type,ion_charge*ion_charge,gkpot_ii,gkR0_ii,gkaR0_ii,R,dEpotent,dQuant);
  // i-e
  return PotGK(type,-ion_charge,gkpot_ei,gkR0_ei,gkaR0_ei,R,dEpotent,dQuant);
}



double PotentialJONES(int type,double R,double& dEpotent, double &dQuant){
  double Rlam=0.;
  switch(type){
  case IONION:
    Rlam=5*R/Lambda_pp;
    break; 
  case ELCELC:
    Rlam=5*R/Lambda_ee;
    break;
  case IONELC:
  case ELCION:
    Rlam=5*R/Lambda_ep;
    break;
  }

  double r12=pow(Rlam,-12);
  double r6=pow(Rlam,-6.);
  dEpotent=E_jones*(r12-r6);
  double k=-E_jones*(-12.*r12+6*r6)/R;
  dQuant=0.;

  if(Pauli_part && (type==ELCELC ||(!non_symm && type==IONION))){ //only for electrons or protons (not for ions)
    double Rlam=R/Lambda_pauli;
    double ex=exp(-Rlam*Rlam);
    dEpotent+=0.5*ex;
    k+=ex*Rlam/Lambda_pauli;
    dQuant+=-0.5*Rlam*Rlam*ex;
  }
  
  dQuant=dEpotent-dQuant;
  return k;
}

 



TableFunction *TABEee=NULL,*TABFee=NULL,*TABQee=NULL;
TableFunction *TABEep=NULL,*TABFep=NULL,*TABQep=NULL;
TableFunction *TABEpp=NULL,*TABFpp=NULL,*TABQpp=NULL;

double PotentialTAB(int type,double R,double& dEpotent,double& dQuant){
  double k=0;
  
  switch(type){
  case IONION:
    dEpotent=(*TABEpp)(R);
    if(TABQpp)dQuant=(*TABQpp)(R);
    k=(*TABFpp)(R);
    break;
  case ELCELC:    
    dEpotent=(*TABEee)(R);
    if(TABQee)dQuant=(*TABQee)(R);
    k=(*TABFee)(R);     
    break;
  case IONELC:
  case ELCION:
    
    dEpotent=(*TABEep)(R);
    dQuant=(*TABQep)(R);
    k=(*TABFep)(R);   
   
    break;
  }

  if(Pauli_part && (type==ELCELC ||(!non_symm && type==IONION))){ //only for electrons or protons (not for ions)
    double Rlam=R/Lambda_pauli;
    double ex=exp(-Rlam*Rlam);
    dEpotent+=0.5*ex;
    k+=ex*Rlam/Lambda_pauli;
    //dQuant+=-0.5*Rlam*Rlam*ex;
  }
  
  return k;
}


void ReadPotential(char *file){
  if(TABEpp && (TABEpp != TABEee) )delete TABEpp;
  if(TABFpp && (TABFpp != TABFee) )delete TABFpp;
  if(TABQpp && (TABQpp != TABQee) )delete TABEpp;

  if(TABEee)delete TABEee;
  if(TABFee)delete TABFee;
  if(TABQee)delete TABQee;
  if(TABEep)delete TABEep;
  if(TABFep)delete TABFep;
  if(TABQep)delete TABEep;
  
  Open_param_file(file);
  if(Set_position("# e-e interaction:")<0)
    serror("Can't find '# e-e interaction:' entry in %s\n",file);
  
  TABEee= new TableFunction(Get_cur_file(),1,2);
  Set_position("# e-e interaction:");
  TABFee= new TableFunction(Get_cur_file(),1,3);
  Set_position("# e-e interaction:");
  TABQee= new TableFunction(Get_cur_file(),1,4);

  if(Set_position("# e-i interaction:")<0)
    serror("Can't find '# e-p interaction:' entry in %s\n",file);
  
  TABEep= new TableFunction(Get_cur_file(),1,2);
  Set_position("# e-i interaction:");
  TABFep= new TableFunction(Get_cur_file(),1,3);
  Set_position("# e-i interaction:");
  TABQep= new TableFunction(Get_cur_file(),1,4);

  int st=Get_stop();
  Set_stop(0);
  if(Set_position("# i-i interaction:")<0){
    TABEpp=TABEee;
    TABFpp=TABFee;
    TABQpp=TABQee;
  }
  else{
    TABEpp= new TableFunction(Get_cur_file(),1,2);
    Set_position("# i-i interaction:");
    TABFpp= new TableFunction(Get_cur_file(),1,3);
    Set_position("# i-i interaction:");
    TABQpp= new TableFunction(Get_cur_file(),1,4);
  }
  Set_stop(st);
  Close_param_file();
}


double PotentialCUT1(int type,double R,double& dEpotent, double &dQuant){
  double R_inv;
  double k;

  switch(type){

  case IONION:
      R_inv=ion_charge/R;
      dEpotent=R_inv*ion_charge;
      k=R_inv*R_inv;
      break;
  case ELCELC:
      R_inv=1./R;
      dEpotent=R_inv;
      k=R_inv*R_inv;
      break;
  case IONELC:
  case ELCION:
    R_inv=1./R;
    if(R>0.3667){
      dEpotent=-R_inv*ion_charge;
      k=R_inv*dEpotent;
    }
    else{
      double y=3.9508*R;
      double sy=sin(y);
      dEpotent=ion_charge*(0.415339*pow(sy,56.222)- 3.);
      k=-92.256*pow(sy,55.222)*cos(y)*ion_charge;
    }
    break;
  }
  dQuant=0.;
  return k;
}







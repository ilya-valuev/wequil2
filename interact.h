# ifndef _INTERACT_H
# define  _INTERACT_H

# include <stdio.h>
# include <string.h>
# include "rparam.h" 
# include "common.h"
 

extern double Lambda_pp,Lambda_ee,Lambda_ep,Lambda_pauli;

extern int Pauli_part, equal_lambda, non_symm;

extern double ion_charge;

extern int neg_cut;
extern double E_cut;
extern double E_negcut;
extern double E_jones;

// Gurski-Krasko parameters (r0, a/r0)
enum GK_POT{
  GK_KELBG, GK_COUL, GK_GK
};
extern double gkR0_ee, gkaR0_ee;
extern double gkR0_ei, gkaR0_ei;
extern double gkR0_ii, gkaR0_ii;
extern int gkpot_ee, gkpot_ei, gkpot_ii;


enum inter_types{
 IONION = 0,
 ELCION = 0x1,
 IONELC = 0x2,
 ELCELC = 0x3
};

typedef double (*EffPotential)(int type,double R,double& dEpotent, double &dQuant);

extern double unit_l;

double PotentialKELBG(int type,double R,double& dEpotent, double &dQuant);
double PotentialERF(int type,double R,double& dEpotent, double &dQuant);
double PotentialCUT(int type,double R,double& dEpotent, double &dQuant);
double PotentialLN(int type,double R,double& dEpotent, double &dQuant);
double PotentialJONES(int type,double R,double& dEpotent, double &dQuant);
double PotentialTAB(int type,double R,double& dEpotent,double& dQuant);
double PotentialCUT1(int type,double R,double& dEpotent, double &dQuant);
double PotentialDEUTSCH(int type,double R,double& dEpotent, double &dQuant);
double PotentialGURKR(int type,double R,double& dEpotent, double &dQuant);
void   WritePotential(const char *file,double start, double L,int n_res,double  f(int,double ,double&,double&));
void ReadPotential(char *file);

double Correction(int type,double potential(int, double, double&,double&),double T);


struct potspec_t{
  EffPotential potential;
  char charpot[51];
  double Clam_ep, Clam_ee;
  double Lambda_set;
  double Lambda;
  int pot_corr;
  double pstart, pend;
  int pn;

  potspec_t():pstart(-1),pend(-1),pn(0){}

  int read_GK_params(int type,int &pot,double &gkaR0,double &gkR0){
    const char *sp= type <0 ? "" : type==0 ? " i-i" : type==0x3 ? " e-e" : type==0x1 ? " e-i" : " i-e";
    char format[200], tmpstr[500];
    sprintf(format,"Gurski-Krasko%s parameters*:> %%s",sp);
    if(Read_param(format,tmpstr)){
      if(strstr(tmpstr,"Kelbg"))
        pot=GK_KELBG;
      else if(strstr(tmpstr,"Coulomb"))
        pot=GK_COUL;
      else if(sscanf(tmpstr,"%lf, %lf",&gkR0,&gkaR0)==2){ // R0, a
        pot=GK_GK;
        gkaR0/=gkR0; // R0, a/R0
      }
      else serror("Unrecognized Gurski-Krasko parameters%s: %s!\n",sp,tmpstr);
      return 1;
    }
    return 0;
  }


  void read_spec(const char *cfgfile){
    int stp=Get_stop();
    char tmpstr[250];
    Set_stop(1);
    Read_param("Potential: %s",tmpstr);
    strncpy(charpot,tmpstr,50);
    if(strstr(tmpstr,"Kelbg"))potential=PotentialKELBG;
    else if(strstr(tmpstr,"Lennard-Johnes"))potential=PotentialJONES;
    else if(strstr(tmpstr,"Gurski-Krasko")){
      potential=PotentialGURKR;
      
      Set_stop(0);
      
       // default
      if(!read_GK_params(3,gkpot_ee,gkaR0_ee,gkR0_ee)) //ee, default is KELBG
        gkpot_ee=GK_KELBG; 
      if(!read_GK_params(0,gkpot_ii,gkaR0_ii,gkR0_ii)) //ii, default is Coulomb
        gkpot_ii=GK_COUL; 
      if(!read_GK_params(0x1,gkpot_ei,gkaR0_ei,gkR0_ei)){ //ei, default is what read from default
        if(!read_GK_params(0x2,gkpot_ei,gkaR0_ei,gkR0_ei)){
          if(!read_GK_params(-1,gkpot_ei,gkaR0_ei,gkR0_ei))
            serror("Can't find any suitable Gurski-Krasko parameters for e-i interaction!\n");
        }
      }
      Set_stop(stp);
      
    }
    else if(strstr(tmpstr,"Deutsch"))potential=PotentialDEUTSCH;
    else if(strstr(tmpstr,"Erf"))potential=PotentialERF;
    else if(strstr(tmpstr,"Cutoff")){
    if(!strcmp(tmpstr,"Cutoff1")){
      potential=PotentialCUT1;
    }
    else{
      potential=PotentialCUT;
      Read_param("Cutoff value*: %lf",&E_cut);
    }
    }
    else if(strstr(tmpstr,"ln")){
      potential=PotentialLN;
      Read_param("Cutoff value*: %lf",&E_cut);
    }
    else if(strstr(tmpstr,"table") && cfgfile){
      potential=PotentialTAB;
      Read_param("Potential table file: %s",tmpstr);
      Close_param_file();
      ReadPotential(tmpstr);
      Open_param_file(cfgfile);
    }
    else serror("Unknown potential type specified!\n");

    Set_stop(0);

    if(Read_param("Pauli part: %s",tmpstr))
      if(strstr(tmpstr,"n"))Pauli_part=0;


    Lambda_set=0.;
    if(Read_param("R0: %s",tmpstr)){
      if(strstr(tmpstr,"default"))Lambda_set=0.;
      else Lambda_set=atof(tmpstr);
    }

    Clam_ep=1.;
    Clam_ee=1.;
    Read_param("e-e R0 coefficient: %lf",&Clam_ee);
    Read_param("e-p R0 coefficient: %lf",&Clam_ep);

    
    
    pot_corr=0;
    if( Read_param("Potential correction: %s",tmpstr)){
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
    }
     
    if(pot_corr)
      printf("The potential will be corrected in %s%s%s parts.\n",
       (pot_corr&0x1 ? "ee" : ""),
       (pot_corr&0x2 ? (pot_corr&0x1 ? ", ii" : "ii"): ""),
       (pot_corr&0x4 ? (pot_corr&0x3 ? ", ei" : "ei"): "") );
    else
       printf("The potential will not be corrected.\n");
  
 
    pn=0;
    if( Read_param("Write potential: %lf, %lf, %d",&pstart, &pend, &pn)==3){
      if(pstart>=pend || pn <2)
        serror("Incorrect potential output settings!\n");
      printf("The potential will be written for the range [%f, %f], %d points\n",pstart,pend,pn);
    }
    else 
      pn=0;
    Set_stop(stp);

  
  }

  void write_pot(const char *pfile, double Lcell){
    WritePotential(pfile,pn? pstart : 0, pn? pend: Lcell/2, pn? pn: 1000,potential);
  }

  void calc_lambda(double T, double ini_Te, double ini_Ti, double mass){
     
    double me=0.9109534e-30;
    double qe=1.602e-19;
    double Kb=1.38e-23;
    double unit_t=qe*qe/(4*M_PI*8.854e-12)*sqrt(me)/pow(Kb*T*1.e4,3./2.);
    double unit_h=Kb*T*1e4*unit_t;
    double h_qwer=1.0545887e-34;

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
      }
      else{
        if(strstr(charpot,"Deutsch")){
          printf("Setting lambda values for Deutsch potential !!!\n");
          Lambda_pauli=h_qwer/(unit_h);
          Lambda=Lambda_pauli/sqrt(2*M_PI);
          double m_ee=0.5, m_pp=0.5*mass, m_ep=mass/(1.+mass);

          Lambda_ee=Lambda/sqrt(m_ee);
          Lambda_pp=Lambda/sqrt(m_pp);
          Lambda_ep=Lambda/sqrt(m_ep);
        }
        else /*if(strstr(Gas.charpot,"Kelbg"))*/{
          printf("Setting lambda values for Kelbg potential !!!\n");

          Lambda_pauli=h_qwer/(unit_h);
          Lambda=Lambda_pauli;
          double m_ee=0.5, m_pp=0.5*mass, m_ep=mass/(1.+mass);
          double t_ep=ini_Te;
          double t_pp=ini_Ti;
          double t_ee=ini_Te;

          Lambda_ee=Lambda/sqrt(2*m_ee*t_ee);
          Lambda_pp=Lambda/sqrt(2*m_pp*t_pp);
          Lambda_ep=Lambda/sqrt(2*m_ep*t_ep);
        }
      }
    }
  }

  void calc_correction(double T){
    // adjusting Gurski-Krasko parameters
    double c =0.53e-10/unit_l; // given in a.u. by the input 
    gkR0_ee*=c;
    gkaR0_ee/=c;  
    gkR0_ei*=c;
    gkaR0_ei/=c;  
    gkR0_ii*=c;
    gkaR0_ii/=c;
    
    if(pot_corr&0x1)
      Correction(ELCELC,potential,T);
    if(pot_corr&0x2)
      Correction(IONION,potential,T);
    if(pot_corr&0x4)
      Correction(IONELC,potential,T);    
  }
};

# endif



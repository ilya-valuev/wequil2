# ifndef ANMFIELD_H
# define ANMFIELD_H

# include "../interact.h"
# include "../analyser.h"
# include "statist.h"

enum {ON_SPACE=0,ON_ION=1,ON_ELC=2};
enum {COULOMB=0x1, EFFECTIVE=0x2, ALL=0x3};

typedef Distribution Distrvect[4];
//typedef Correlation  Corrvect[3];

class AnMField: public Analysator{
  int mc_np, mc_rnd, mc_nop;
  Vector_3 *mc_tpx;
  char **mc_tpname;
  double mc_r0, mc_r1;
  int mc_nr, mc_xmgr, mc_prcoord, mc_prav,mc_correl;
  cList* mc_plist;
  double mc_F0;
  int mc_type, mc_force, mc_ntype /* ion*/;
  double mc_cut;
  int mc_from;
  int mc_anim;
  int mc_grid;

  Distrvect *Mf;
  Distrvect *Ef;
  Correlation *MfCorr;
  Correlation *EfCorr;
  //Distrvect Mfav;
  
  EffPotential mc_potential;

  void output_mfdistr(char *filename, char *dname, Distrvect *Mf);
  //int output_vcorr(char *fileform,Correlation* Corr, int npart);
  void output_mfcorr(char *filename,Correlation* MfCorr);
public:

  AnMField(): mc_tpx(NULL), Mf(NULL), Ef(NULL), MfCorr(NULL), EfCorr(NULL){
    mc_np=0;
    mc_rnd=0;
    mc_nop=0;
    mc_xmgr=0;
    mc_prcoord=0;
    mc_correl=0;
    mc_force=0x1;
    mc_ntype=0 /* ion*/;
    mc_plist=NULL;
    mc_tpname=NULL;
    mc_cut=-1;
    mc_from=0;
    mc_anim=0;
    mc_grid=-1;
  };
  ~AnMField(){
    if(mc_tpx && mc_type==ON_SPACE)
      delete [] mc_tpx;
    if(Mf)delete [] Mf;
    if(Ef)delete [] Ef;
    if(MfCorr)delete [] MfCorr;
    if(EfCorr)delete [] EfCorr;
    if(mc_plist)delete mc_plist;
    
    if(mc_tpname){
      int i;
      for(i=0;i<mc_np;i++){
        delete [] mc_tpname[i];
      }
      delete [] mc_tpname;
    }
  }
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  //int ProcessSplit();
  int Process(char *name);  
};

# endif

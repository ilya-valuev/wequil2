# ifndef ANVCORR_H
# define ANVCORR_H

# include "../analyser.h"
# include "statist.h"



class AnVcorr: public Analysator{
  Correlation *vc_Ci, *vc_Ce;
  int vc_ion, vc_elc;
  int vc_type;
  int vc_ne,vc_ni, vc_sav;
  int nesteps, nisteps;

  int output_vcorr(char *fileform,Correlation* Corr, int npart);
public:

  AnVcorr(){
    vc_ion=0;
    vc_elc=0;
    vc_type=VC_FFT;
  };
  ~AnVcorr();
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);  
};



# endif








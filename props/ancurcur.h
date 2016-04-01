# ifndef ANCURCUR_H
# define ANCURCUR_H

# include "../analyser.h"
# include "statist.h"


class AnCurCur: public Analysator{
  Correlation *vc_Ci, *vc_Ce, *vc_Ct, *vc_Cei, *vc_Cie;
  int vc_ion, vc_elc, vc_tot, vc_cross;
  int vc_type;
  int vc_ne,vc_ni, vc_sav;
  int vc_norm;
  int nesteps, nisteps;

  int output_vcorr(char *fileform,Correlation* Corr);
  void allocate_data();
  void init_data();
  void prepare_total();
  void insert_end();
  void do_correlations();
public:
  AnCurCur(){
    vc_ion=0;
    vc_elc=0;
    vc_tot=0;
    vc_cross=0;
    vc_norm=0;
    vc_type=VC_FFT;
  };

  ~AnCurCur();
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);
};



# endif








# ifndef ANDSTRUC_H
# define ANDSTRUC_H

# include "../analyser.h"
# include "../statist.h"




typedef Correlation *CorrP;


# define NCORR 5  // 0-i, 1-e, 2-ie, 3-ei, 4-sum (for check only)

extern char *AnDsymb[];



enum { FOUT_NONE=0, FOUT_R=0x1, FOUT_I=0x2, FOUT_M=0x4 };


class AnDstruc :public Analysator{

  int nd_n, nd_ni, nd_type, nd_nvect;
  Vector_3 *nd_kvect;
  //cList *klist;
  //int *
  CorrP nd_r[NCORR];
  float factor_i, factor_e;
 
  int nesteps, nisteps, nsteps;
  float cur_isr, cur_isi;

  int wr_tc, wr_fc;

  int f_corr_int;
  float f_ws, f_we, f_dw;
  float f_itime;
  int f_out[6];
  int f_outn;
  // smoothing window for Fourrier transform, in Wpe
  float out_w;

  int fileout(char *fileform);
public:

  AnDstruc(){
    nd_nvect=0;
    nd_type=VC_FFT;
  };

  ~AnDstruc();
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);  
};



# endif








# ifndef ANTHERM_H
# define ANTHERM_H


# include <stdio.h>
# include "../analyser.h"
# include "statist.h"



class AnTherm: public Analysator{
  double Te, Ti, T;
  RunAverage STe, STi, ST;
  Statistics ATe, ATi, AT;
  double tav, tavi;
  double t0, t0e, t0i, t0iw;
  FILE *f, *fi;
  int one_center;
  int draw_av;
public:
  AnTherm(): STe(&Te),STi(&Ti), ST(&T) {}
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int Process(char *name);  
};


# endif

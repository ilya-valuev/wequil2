# ifndef ANTHERM_H
# define ANTHERM_H


# include <stdio.h>
# include "../analyser.h"
# include "../statist.h"



class AnForce: public Analysator{
  Measurement m[20];
  /*double values[10];
  double rav[10];
  double &MaxR, &MinR, &MaxF, &MinF, &MaxV, &MinV, &Uii, &Uei, &Uee, Te, Ti, T;
  RunAverage STe, STi, ST; */
  double tav, tavi;
  double t0, t0e, t0i, t0iw;
  FILE *f, *fi;
public:
  AnForce(){}
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int Process(char *name);
};


# endif

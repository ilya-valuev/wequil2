# ifndef ANPAIR_H
# define ANPAIR_H
       
# include "../analyser.h"
# include "statist.h"
 
class AnPair: public Analysator{
  double rr_coeff;
  Distribution DRRee,DRRep,DRRpp, DRReep, DRReeu;
  int spin;
  int WriteDRR(char *file);
public:
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);  
};


# endif

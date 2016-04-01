# ifndef ANVDISTR_H
# define ANVDISTR_H

# include "../analyser.h"
# include "statist.h"

// how to calculate v
# define BASE_ABS       0    // as is
# define BASE_COMMON    1    // substract common center-of-mass velocity
# define BASE_SEPARATE  2    // substract center of mass velocity for each component


class AnVdistr: public Analysator{
  Distribution DVe[4],DVp[4];
  double Te, Ti;
  Statistics STe, STi;
  int max_comp;
  double vk;
  int write_file(char *fname, Distribution *DV,double mass, double T);
  int av_coord;
  int distr_base;
public:

  AnVdistr(): STe(&Te), STi(&Ti){};
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);  

};


double maxwell_vx(double m, double kT,double v);

double maxwell_v2(double m, double kT,double v);


# endif

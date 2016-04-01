# ifndef ANTRJLOG_H
# define ANTRJLOG_H
       
# include "../analyser.h"
# include "statist.h"
 
class AnTrjLog: public Analysator{
  FILE *fp;
  int init_file();
  char name[2][100];
  double scale;
public:
  int Init(AnalyseShell *parent,int initype=0);
  int Step();
  int ProcessSplit();
  int Process(char *name);  
};


# endif

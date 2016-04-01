# ifndef ANALYSE__H
# define ANALYSE__H

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include "vector_3.h"
#include "common.h"
#include "rparam.h"
#include "statist.h"

# include "plasma.h"
/*
enum dtypes {
 NOTHING =0,
 COORD   =0x1,
 VEL     =0x2,
 FLOW    =0x4 };*/

int Init_analyse(char *cfgfile,int nn,double L, double Gamma, double T,char *daname,double tstep);
void Step_analyse(Vector_3 *xx, Vector_3 *v);
void Process_analyse();

int InitTrajectory(long &type,long step,char *datafile,char *dname,double &anL,double &anGamma,double &anT, double &animass);

void Refine_step();
long gettypesize(long type, int n);

void WriteDRR(char *file);


void Flow(int n, Vector_3 *v, Vector_3* flowe, Vector_3* flowp);

extern Distribution DRRee,DRRep,DRRpp;


extern FILE *tFile, *new_trj;
extern double dt;
extern long nsteps,ref_type;
extern int npart;
extern double animass;
extern Vector_3 *anx, *anv;


#endif








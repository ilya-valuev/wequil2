#ifndef _PCYCLE_H
# define _PCYCLE_H




typedef struct {
 double p1;
 double p2;
 int i;
 int n;
 int spec;
 char s[16];
} Parameter;


# define DEFAULT  0
# define GNU      1

extern int tblform;


void SetTableform(int tf);


double get_pc(Parameter& p);


int InitParameters(Parameter** pp);

int CheckData(char *logfile,char *datafile,char *dataname,int np,Parameter* p, int ask=1);

long CycleCount(int np,Parameter* p);

long CurrentCount(int np, Parameter *p);


void BeginFrame(FILE *f1,char *dataname,int np,Parameter *p);


void MiddleFrame(FILE *f1,int np,Parameter *p);


void WriteStatus(int done,char *sfile,char *datafile,int np,Parameter *p);


void StatusLine(char *str,int np,Parameter *p);


int StopStatus(char *file, int stat);

void StartTime(void);

# endif

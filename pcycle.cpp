# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <fcntl.h>

# ifndef UNIX
# include <io.h>
# endif

#include <time.h>
#include <sys/types.h>


# include "pcycle.h"
# include "rparam.h"
#include "common.h"


int tblform=DEFAULT;

void SetTableform(int tf){
  tblform=tf;
}

time_t prg_time=0;
clock_t prg_clock=0;

void StartTime(void){
 time(&prg_time);
 prg_clock=clock();
}

long GetTime(void){
 time_t c_time;
 time(&c_time);
 return (long)(c_time-prg_time);
}


float GetClock(void){
 clock_t c_time;
 c_time=clock();
 return ((float)(c_time-prg_clock))/CLOCKS_PER_SEC;
}


double get_pc(Parameter& p){
 if(p.n>1)return p.p1+(p.p2-p.p1)*p.i/(p.n-1);
 else return p.p1;
}

int InitParameters(Parameter** pp){
 int np,i;
 char str[256];
 Read_param("Number of parameters: %d",&np);

 Parameter *p;
 *pp=p=new Parameter[np];
 if(!p)serror("Error allocating memory for parameters!\n");
 printf("Parameter: n=%d\nValues:\n",np);
 for(i=0;i<np;i++){
  sprintf(str,"p%d*: %%lf,%%lf,%%d",i+1);
  int c=Read_param(str,&p[i].p1,&p[i].p2,&p[i].n);
  if(c==0 || c==2)serror("Can't read parameter #%d value(s)!\n",i+1);
  if(c<2){
   p[i].n=1;
   p[i].p2=p[i].p1;
  }
  strncpy(p[i].s,AcFormat,15);
  p[i].s[15]=0;
  p[i].i=0;
  printf("%f, %f, %d\n",p[i].p1,p[i].p2,p[i].n);
 }
 return np;
}


int careful=1;

int CheckData(char *logfile,char *datafile,char *dataname,int np,Parameter* p, int ask){
 int no_remove=0,i;
 FILE *f1;
 char str[256];
 f1=fopen(logfile,"rt");
 if(f1){
  fclose(f1);
  Open_param_file(logfile);
  Set_stop(0);
  Read_param("Data: %s",str);
  Set_stop(1);
  if(!strcmp(str,dataname)){ // continuation
   Read_param("Finished: %s",str);
   if(!strcmp(str,"no")){
    no_remove=1;
    for(i=0;i<np;i++){
     sprintf(str,"p%d : %%d",i+1);
     if(Read_param(str,&p[i].i)==0)serror("Can't read parameter #%d cycle count!\n",i+1);
    }
    printf("\n %s calculation continued from ",datafile);
    for(i=0;i<np;i++){
     if(i!=0)printf(",");
     printf("%d",p[i].i);
    }
    printf(".\n");
   }
  }
  Close_param_file();
 }

 
 if((f1=fopen(datafile,"rt"))!=NULL){
   fclose(f1);
   if(!no_remove){     
     if(ask){
       printf("File %s already exists. Overwrite? (y/n)\n",datafile);
       if(getchar()=='n')exit(1);
       remove(datafile);
     }
     if(careful){
        serror("File %s already exists with 'finished' or undefined status.\n"
               "Please remove it or use -new option to proceed.\n",datafile);
     }
     else remove(datafile);
   }
 } 
 else if(no_remove)no_remove=2;

 return no_remove;
}


long CycleCount(int np,Parameter* p){
  long lm=1,sum=0;
  int i;
  (p[0].i)++;
  for(i=0;i<np;i++){
   if(p[i].i>=p[i].n && i!=np-1){
    p[i].i=0;
    (p[i+1].i)++;
   }
   sum+=lm*p[i].i;
   lm*=p[i].n;
  }
  if(p[np-1].i<p[np-1].n)return sum+1;
  else return 0;
}


long CurrentCount(int np, Parameter *p){
  long lm=1,sum=0;
  int i;
  
  for(i=0;i<np;i++){
   sum+=lm*p[i].i;
   lm*=p[i].n;
  }
  return sum+1;
} 


void BeginFrame(FILE *f1,char *dataname,int np,Parameter *p){
 int i;
 char str[256];
 char delimet[256]="";
 if(tblform==GNU)strcpy(delimet,"# ");

 fprintf(f1,"\n\n%sData name: %s\n%sParameters:\n",delimet,dataname,delimet);
 for(i=0;i<np;i++){
  if(p[i].n!=0)fprintf(f1,"%s%s : %f, %f, %d\n",delimet,p[i].s,p[i].p1,p[i].p2,p[i].n);
  else fprintf(f1,"%s%s : %f\n",delimet,p[i].s,p[i].p1);
  fprintf(f1,"%sValues: ",delimet);
  int c=p[i].i;
  if(p[i].n!=0){
   for(;p[i].i<p[i].n;(p[i].i)++)fprintf(f1,"(%d,%f);",p[i].i,get_pc(p[i]));
  }
  else{
   fprintf(f1,"(%d,%f);",0,p[i].p1);
   p[i].n=1;
  }
  p[i].i=c;

  fprintf(f1,"\n");
 }

 sprintf(str,"\n%%%d.%ds",5*np+3,5*np+3);
 fprintf(f1,str,strcat(delimet,"Cycle counts "));

 for(i=0;i<np;i++){
  fprintf(f1,"%12.12s ",p[i].s);
 }
}

void MiddleFrame(FILE *f1,int np,Parameter *p){
  int i;
  char str[256];

  if(tblform==GNU)fprintf(f1," ");
  else fprintf(f1,"(");

  for(i=0;i<np;i++){
   fprintf(f1,"%4d",p[i].i);
   if(i!=np-1){
     if(tblform==GNU)fprintf(f1," ");
     else fprintf(f1,",");
   }
  }

  if(tblform==GNU)fprintf(f1,"   ");
  else fprintf(f1,"): ");

  str[0]=0;
  for(i=0;i<np;i++){
   sprintf(str+strlen(str),"%12.4f ",get_pc(p[i]));

  }
  fprintf(f1,"%s",str);
}



# ifdef UNIX
# include <unistd.h>
# endif



void WriteStatus(int done,char *sfile,char *datafile,int np,Parameter *p){
  FILE *f1;
  int i;
  f1=Err_fopen(sfile,"wt");
  fprintf(f1,"Current status file\n");
# ifdef UNIX
  char str[200];
  gethostname(str,200);
  long pid=(long)getpid();
  fprintf(f1,"Host: %20s\nPID: %ld\n",str,pid);
# endif
  fprintf(f1,"Elapsed time: %lds\n",GetTime());
  fprintf(f1,"Processor time: %fs\n",GetClock());
  fprintf(f1,"Data: %s\n\n",datafile);



  if(!done){
    fprintf(f1,"Finished: no\n");
    fprintf(f1,"Cycle counts:\n");
    for(i=0;i<np;i++){
      fprintf(f1,"p%d : %d\n",i+1,p[i].i);
    }
  }
  else fprintf(f1,"Finished: yes\n");
  fprintf(f1,"\n");
  fclose(f1);
}

void StatusLine(char *str,int np,Parameter *p){
 sprintf(str,"Cycle: ");
 int i;
 for(i=0;i<np;i++){
   sprintf(str+strlen(str),"%d/%d ",p[i].i+1,p[i].n);
 }
}

int StopStatus(char *file, int stat){
  FILE *f1;
  f1=_fdopen(_open(file,O_RDWR),"r+");
  if(!f1)serror("Can't reopen file %s",file);

  char str[250];
  fseek(f1,0,SEEK_SET);
  //printf("Pos: %ld\n",ftell(f1));
  fscanf(f1,"%s",str);
  int stp=0;
  if(strstr(str,"stop"))stp=1;
  //printf("Stop stat: %s\n",str);

  fseek(f1,0,SEEK_END);

  long i=ftell(f1);
  int k=0;
  int num, di;
# ifdef UNIX
  num=2;
  di=1;
# else
  num=4;
  di=2;
# endif

  do{

    i--;

    fseek(f1,i,SEEK_SET);
    if(getc(f1)=='\n')k++;
  }while(k<num && i>0);

  fseek(f1,i+di,SEEK_SET);
  fprintf(f1,"Proceed: %d  %lds\n",stat,GetTime());

  fclose(f1);
  return stp;
}

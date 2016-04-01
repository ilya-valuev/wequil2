# ifndef PLASMAR_H
# define PLASMAR_H

# include "record.h"
# include "plasma.h"


class PlasmaRecord : public Record{
  Plasma *plasmap;
  int wne;
  int wni;
  int wq;
  double wT;
  double wGamma;
  double wstp;
  double wmass;
  char wpot[50];
  int wnsymm;

  

  long wenseq,wions, welec;

  double dstep;
  double real_step;
protected:
  char wdname[100];
  void write_charhead();
 
public:
  int valid;
  int wtype;
  PlasmaRecord(){
    valid=0;
    wtype=0;

    dstep=1.;
    real_step=0.;
    plasmap=NULL;
  }
  int Open(char *filename);
  int Check(char *file,int type,Plasma &Gas,long wr_int,
	    long wr_ions,long wr_enseq,int num_only=0);
  int Init(char *file,int type,Plasma &Gas,long wr_int,
	    long wr_ions,long wr_enseq);

  double AdjustInterval(double dt);
  int Query();
  int Truncate(long new_nstp){
    if(new_nstp<1 || new_nstp>nsteps)return 0; // invalid 

    sleep_FP();
    curpos=step_position(new_nstp);
    wake_FP();
    write_stepmark(new_nstp); // stepmark specifies unfinished step
    sleep_FP();
    return 1;
  }
  long ReloadGas(Plasma &Gas, int adjfile=0);
    
  long Step(){ // step converter
    real_step+=dstep;
    long stp=(long)(real_step);
    //printf("stp: %ld, real=%f\n",stp,real_step);
   
    for(;step+1<=stp;){
      //printf("rstp: %ld\n",step);
      Record::Step();
      if(state==R_WAITQ){ // control correction
	real_step=((double)step+1e-5)/dstep;
	return step;
      }
    }
    return step;
  }


  double file_dt(){
    return wstp;
  }

  double file_edt(){
    return welec*wstp;
  }

  double file_idt(){
    return wions*wstp;
  }

  double file_edelay(){
    if(wenseq<0)return 0.;
    else return wenseq*wstp;
  }

  int AllocPlasma();

  int GetAllNext();

  friend class Analysator;
  friend class AnalyseShell;
};

void WriteTypeOut(char *str,int wtype);
int WriteTypeIn(char *str);
    



// data status specificators
enum{
  PR_UNKNOWN=0,
  PR_ELCVEL=1,
  PR_ELCCOORDS=2,
  PR_IONVEL=4,
  PR_IONCOORDS=8,
  PR_FLOW=16
};


    

#endif







# ifndef ANALYSER_H

# define ANALYSER_H

# include <string.h>
# include "plasmar.h"


double random2(void);



// types of time split
enum {
  NO_SPLIT =0,
  AV_SPLIT =0x1,
  PROF_SPLIT=0x2
};


class Analysator;

/*
int ANelements=2;
AnPair
Analysator *AnArray[]={*/
  
  
# define MAX_ANALYSERS  20


class AnalyseShell{
  char out_dir[256];
protected:
  
  double dt;
  double edt, idt, edelay;
  double max_time;
  double start_time;
  int trjstep;

  char dname[100];

  PlasmaRecord *rec;
  Vector_3 *anx, *anv;
  int xx_used;
  Plasma *gasp;
  int ivel, icoord, evel, ecoord;

  Analysator *a[MAX_ANALYSERS];
  int na;

  int sp_type;

  double sp_time;

  double sp_time_act;
  double sp_time_act_next;

  double t0s;
  int sp_num;

  long sp_nsteps;
  long sp_nsteps_next;

public:
  int lastsplit;
  double rtime;

  long nsteps;
  long cur_rp; // record pos
  long start_rp; // starting record pos
  int status; // current status
  int trj_end;  // end indicator

  int update;

  int valid;

  AnalyseShell(){
    rec=NULL;
    gasp=NULL;
    nsteps=1;
  }

  ~AnalyseShell();
  

 int Setrec(PlasmaRecord *pr){
    rec=pr;
    if(!rec->valid)return 0;
    if(!rec->AllocPlasma())return 0;
    gasp=pr->plasmap;
    anx=gasp->x;
    anv=gasp->v;
    xx_used=0;
   
    dt=rec->file_dt()*gasp->Wpe;
    edt=rec->file_edt()*gasp->Wpe;
    idt=rec->file_idt()*gasp->Wpe;
    edelay=rec->file_edelay()*gasp->Wpe;

    printf("dt: %f\n",dt);
    nsteps=rec->getNSteps();     

    ivel=rec->FrameWith("ion_veloc");
    icoord=rec->FrameWith("ion_coord");
    evel=rec->FrameWith("electron_veloc");
    ecoord=rec->FrameWith("electron_coord");
  
    if((rec->wdname)[0]!=0)strcpy(dname,rec->wdname);

    na=0;
    valid=1;
    return 1;
  }
  void SetOutDir(char *dir);

  int Init(char *file);


  int AddAnalysator(Analysator *a);

  int StartRead(){
    return rec->FixReadStep(cur_rp*trjstep);
  }
  
  int StepRec();
  int Process();

  double ndone(){
    return ((double)(cur_rp-(start_time)/dt))/nsteps;
  }

  long regsteps(int fr_type);

  long regsteps_next(int fr_type);

  double split_length(){
    if(sp_type!=NO_SPLIT){
      return sp_time_act;
    }
    return max_time-start_time;
  }  

  double split_length_next(){
    if(sp_type!=NO_SPLIT){
      return sp_time_act_next;
    }
    return max_time-start_time;
  }  

  char *GetDataName(){
     return dname;
  }
  friend class Analysator;
  friend class AnPair;
  friend class AnTherm;
  friend class AnVdistr;
  friend class AnVcorr;
  friend class AnDstruc;
  friend class AnCurCur;
  friend class AnMField;
  friend class AnTrjLog;
};


enum{
  COLP_SUM=0,
  COLP_STAT=1 
};  


class col_desc{
public:
  double *col;
  double *nav;
 
  int colprop;

  col_desc(){
    col=NULL;
    nav=NULL;
    
  }

  int init(int np);

  ~col_desc(){
    if(col){
      delete [] col;
      delete [] nav;
      
    }
  }
    
};


class OutputArr{
  int UpdateColumnIR(int cn, TableFunction *data, TableFunction *xn=NULL);
  int UpdateColumnIR(int cn,double F(double));
  double arg_coeff;
public:
  char file[250];
  col_desc *c;
  int ncol;
  int np;
  double ts,te,dt;
  int ftype;
  Pool pxx;
  Pool **pyy, **pnn;

  OutputArr(double t1, double t2, int npoints, int nc):pxx(NULL,sizeof(double)),arg_coeff(1.){
    c=NULL;
    init(t1,t2,npoints,nc);
  }

  OutputArr(int nc):pxx(NULL,sizeof(double)),arg_coeff(1.){
    c=NULL;
    init(nc);
  } 

  int init(double t1, double t2, int npoints, int nc);

  int init(int nc);

  int NextPoint(double x, double *y, double *n=NULL);

  void SetArgCoeff(double c=1.){
    arg_coeff=c;
  }

  void dealloc();
		

  ~OutputArr(){
    dealloc();
  }

  int UpdateColumn(int cn, TableFunction *data, TableFunction *xn=NULL);
  int UpdateColumn(int cn,double F(double));

  int PrintToFile(FILE *f);

  int Load();

};




class Analysator{
protected:
  Pool pOut;
  OutputArr *Out;

  AnalyseShell *psh;
  int valid;
public:
  Analysator(): pOut((void **)&Out,sizeof(OutputArr)) {
    psh=NULL;
    valid=0;
  }

  virtual int Init(AnalyseShell *parent,int initype=0){
    psh=parent;
    valid=0;
    if(!psh || !psh->valid)return 0;
    valid=1;
    return 1;
  }
  virtual int Step(){
    return valid;
  }
  virtual int Process(char *name){
    valid=0;
    return 1;
  }

  virtual int ProcessSplit(){
    return 1;
  }

  virtual ~Analysator(){
    //printf("oooe\n");
  }


  int ExistsOutput(int o_num){
    if(o_num>=0 && pOut.ind > o_num)return 1;
    return 0;
  }

  int AddOutput(double t1, double t2,int npoints, int nc){
    OutputArr *tmp=new OutputArr(t1,t2,npoints,nc);
    pOut.add(tmp);
    // Out[pOut.ind-1].init(t1,t2,npoints,nc);
    return pOut.ind;
  }

  /* irregular output*/
  int AddIrOutput(int nc){
    OutputArr *tmp=new OutputArr(nc);
    pOut.add(tmp);
    // Out[pOut.ind-1].init(t1,t2,npoints,nc);
    return pOut.ind;
  }


  int AssociateOutput(char *filename, int o_num, int reload=0){
    if(!ExistsOutput(o_num))return -1;

    strcpy(Out[o_num].file,filename);

    if(reload){
      if(!Out[o_num].Load()){
	eprintf("Using new '%s'.\n",filename);
      }
    }
    return 1;
  }

  int UpdateFile(int o_num, char *header=NULL){
    if(!ExistsOutput(o_num))return -1;
    FILE *f;
    f=fopen(Out[o_num].file,"wt");
    //printf("ohohoho: %s\n",header);
    if(!f){
      msg_error("Analysator:"
		"UpdateFile:Can't open file '%s'\n",Out[o_num].file);
      return -1;
    }
    if(header){
      fprintf(f,"%c%s\n",'#',header);
    }

    Out[o_num].PrintToFile(f);
    fclose(f);
    return 1;
  }



};

class Measurement:public Analysator{
protected:
  char name[250];
  int dim; // dimension
public:
};

class Evolution: public Measurement{
protected:
  RunAverage *rav;
  int w; // window
public:
  //double unit;
  int GetWindow() const { return w; }
  int SetWindow(int win){
    int i, res=0;
    for(i=0;i<dim;i++)res+=rav[i].set_window(w);
    return res;
  }
  double value(int i=0){
    return rav[i].rav();
  }

};





# endif



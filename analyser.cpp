# include <stdio.h>
# include <ctype.h>

#include <errno.h>
# ifdef UNIX
#include <sys/stat.h>
#include <sys/types.h>
# else
# include <direct.h>
# endif



# include "plasma.h"
# include "rparam.h"
# include "common.h"
# include "analyser.h"

# include "props/anpair.h"
# include "props/antherm.h"
# include "props/anvdistr.h"
# include "props/anvcorr.h"
# include "props/andstruc.h"
# include "props/ancurcur.h"
# include "props/anmfield.h"
# include "props/antrjlog.h"

double random2(void){
  return 1.-2.*((double)rand())/RAND_MAX;
}


void AnalyseShell::SetOutDir(char *tmpstr){
  strncpy(out_dir,tmpstr,255);

  int n=strlen(out_dir);
  if(n>0){
    if(out_dir[n-1]!='/')strcat(out_dir,"/");
  }
}


int AnalyseShell::Init(char *file){
  char tmpstr[256];

  if(!gasp)return 0;

  Set_comment_char(';');
  Open_param_file(file);
  Set_stop(0);

//# ifdef UNIX

  if(Read_param("Output directory: %s",tmpstr)){
    SetOutDir(tmpstr);
# ifdef UNIX
    if(mkdir(out_dir,S_IRWXU|S_IRGRP|S_IROTH)){
# else
    if(_mkdir(out_dir)){
# endif
      if(errno!=EEXIST){
	      msg_error("Can not create directory: %s.\n%s\n",out_dir,strerror(errno));
	      return 0;
      }
    }
  }
//# else
//  SetOutDir("");
//# endif

  max_time=-1;
  start_time=0;
  trjstep=1;

  
  if(!Read_param("Stop  time: %lf",&max_time)){
    if(Read_param("Stop  step: %lf",&max_time))
      max_time*=dt;
  }
  if(!Read_param("Start time: %lf",&start_time)){
    if(Read_param("Start step: %lf",&start_time))
      start_time*=dt;
  }

  //max_time/=gasp->Wpe;
  //start_time/=gasp->Wpe;

  Read_param("Trajectory step: %d",&trjstep);


  if(max_time<0 || max_time > dt*nsteps )max_time=dt*nsteps;
  if(start_time<0)start_time=0;

  if(max_time<start_time){
    msg_error("Start and/or stop time are not consistent"
	   " with trajectory length.\n");
    return 0;
  }

  nsteps=(long)((max_time-start_time)/dt);
  nsteps/=trjstep;
  dt*=trjstep;

  if(!Read_param("Split time: %lf",&sp_time)){
    sp_time=max_time-start_time;
    sp_type=NO_SPLIT;
  }
  else{
    sp_type=NO_SPLIT;
    if(Read_param("Split type:>",tmpstr)){
      if(strstr(tmpstr,"profile"))sp_type|=PROF_SPLIT;
      if(strstr(tmpstr,"ave"))sp_type|=AV_SPLIT;
    }
    else{
      printf("Analyse: Init: cannot read split type, using profile.\n");
      sp_type=PROF_SPLIT;
    }
  }
   
  update=0;
  if(Read_param("Update average data: %s",tmpstr)){
    if(strstr(tmpstr,"y"))update=1;
  }

  

  strcpy(dname,"a");
  Read_param("Set data name: %s",dname);

  cur_rp=(long)((start_time)/dt);
  start_rp=cur_rp;
  trj_end=0;
  rtime=cur_rp*dt;

  t0s=rtime;
  sp_num=0;

  sp_time_act=fmin(sp_time,max_time-start_time);
  sp_time_act_next=sp_time;

  sp_nsteps=(long)(sp_time_act/dt);
  sp_nsteps_next=sp_nsteps;
  
  lastsplit=0;

  Analysator *antmp;
  na=0;
 
  Set_stop(0);
  if(Read_param("Pair distribution function: %s",tmpstr)){
    if(tmpstr[0]=='y'){
      antmp=new AnPair();
      if(!antmp){
	msg_error("Error creating analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
    }
  }

  Set_stop(0);
  if(Read_param("Temperature log: %s",tmpstr)){
    if(tmpstr[0]=='y'){
      antmp=new AnTherm();
      if(!antmp){
	msg_error("Error creating analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
    }
  }

  Set_stop(0);
  if(Read_param("Velocity distribution function: %s",tmpstr)){
    if(tmpstr[0]=='y'){
      antmp=new AnVdistr();
      if(!antmp){
	msg_error("Error creating analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
    }
  }

  Set_stop(0);
  if(Read_param("Velocity autocorrelation: %s",tmpstr)){
    if(tmpstr[0]=='y'){
      antmp=new AnVcorr();
      if(!antmp){
	msg_error("Error creating analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
    }
  }

  Set_stop(0);
  if(Read_param("Dynamic structure factor: %s",tmpstr)){
    if(tmpstr[0]=='y'){
      antmp= new  AnDstruc();
      if(!antmp){
	msg_error("Error creating analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
    }
  }

  Set_stop(0);
  if(Read_param("Current-current corr*: %s",tmpstr)){
     if(tmpstr[0]=='y'){
      antmp= new  AnCurCur();
      if(!antmp){
	msg_error("Error creating current-current analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
     }
  }

  Set_stop(0);
  if(Read_param("Microfield analysis: %s",tmpstr)){
     if(tmpstr[0]=='y'){
      antmp= new  AnMField();
      if(!antmp){
	msg_error("Error creating microfield analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
     }
  }

  Set_stop(0);
  if(Read_param("Animation: %s",tmpstr)){
     if(tmpstr[0]=='y'){
      antmp= new  AnTrjLog();
      if(!antmp){
	msg_error("Error creating animation analyser!\n");
	return 0;
      }
      AddAnalysator(antmp);
     }
  }


  Close_param_file();

  return 1;
}



int AnalyseShell::AddAnalysator(Analysator *an){ 
  if(na>=MAX_ANALYSERS-1){
    msg_error("Too many analysers in one run!\n");
    return 0;
  }

  a[na]=an;
  if(a[na]->Init(this,0)==0){
    msg_error("Error by analyser initialization (#%d)!\n",na+1);
    return 0;
  }
  na++;
  return 1;
}


AnalyseShell::~AnalyseShell(){
  int i;
  for(i=0;i<na;i++){
    //printf("del %d\n",i);

    delete a[i];
  }
}


int AnalyseShell::StepRec(){
  long stp=cur_rp*trjstep;

  int i;


  status=0;
  if(trj_end || (cur_rp-start_rp) > nsteps || stp>=rec->getNSteps()){
    trj_end=1;
    return 0;
  }

# if 0
  status=rec->GetAllNext();
# endif


# if 1
  //printf("stp%ld (%ld) %d ",stp,nsteps,na);
  // reading Gas
  if(rec->Get(stp,ivel,anv)){
    status|=PR_IONVEL;
    //printf("ivel\n");
  }
  if(rec->Get(stp,evel,anv+gasp->ni)){
    status|=PR_ELCVEL;
    //printf("evel");
  }
  if(rec->Get(stp,icoord,anx)){
    status|=PR_IONCOORDS;
    //printf("icoord\n");
  }
  if(rec->Get(stp,ecoord,anx+gasp->ni)){
    status|=PR_ELCCOORDS;
    //printf("ecoord");
  }
# endif

  xx_used=0;

  rtime=cur_rp*dt;
  cur_rp++;


  for(i=0;i<na;i++){
    //printf("dodo %d\n",i );
    a[i]->Step();   // performing analyser steps
  }

  if(rtime-t0s>sp_time_act){ //new split begins
    printf("sp: %d\n",sp_num);

    int lspl;
    if(max_time<rtime+sp_time){
      lspl=1;
      sp_time_act_next=max_time-rtime;
    }
    else{
      lspl=0;
      sp_time_act_next=sp_time;
    }


    // nsteps for the next split
    sp_nsteps_next=(long)(sp_time_act_next/dt);

    
    for(i=0;i<na;i++){
      a[i]->ProcessSplit();
    }

    sp_nsteps=sp_nsteps_next;
    sp_time_act=sp_time_act_next;
    t0s=rtime;
    sp_num++;
    lastsplit=lspl;
  }

  for(i=0;i<trjstep;i++){
    rec->NextReadStep();
  }

  //printf("ee -%d\n",trj_end);
  return status;
}


int AnalyseShell::Process(){
  status=0;
  int i;
  for(i=0;i<na;i++){
    a[i]->Step();   // performing last analyser steps
    lastsplit=1;
    a[i]->ProcessSplit(); // ending splits
  }
  

  for(i=0;i<na;i++){
    a[i]->Process(dname);   // performing last analyser processing
  }
  return 1;
}


long AnalyseShell::regsteps(int fr_type){
  long res=-1;
  long frame_fq=1;
  if(rec){
    int frm=rec->FrameWithSpec(fr_type);
    if(frm<0)return 0;

    res=rec->NRegSteps(cur_rp*trjstep,frm);
    frame_fq=rec->FrameFreq(frm);
  }
  long sp_stp=sp_nsteps/frame_fq;


  if(res>0 && sp_type!=NO_SPLIT){
    if(res<sp_stp)return res;
  }
  return sp_stp;
}
    

long AnalyseShell::regsteps_next(int fr_type){
  long res=-1;
  long frame_fq=1;
  if(rec){
    int frm=rec->FrameWithSpec(fr_type);
    if(frm<0)return 0;
    
    res=rec->NRegSteps(cur_rp*trjstep,frm);
    frame_fq=rec->FrameFreq(frm);
  }
  long sp_stp=sp_nsteps_next/frame_fq; 


  if(res>0 && sp_type!=NO_SPLIT){
    if(res<sp_stp)return res;
  }
  return sp_stp;
}




int col_desc::init(int np){
  col=new double[np];
  nav=new double[np];
  

  if(!col || ! nav){
    msg_error("col_desc: OutputArr: MAE\n");
    return 0;
  }
  colprop= COLP_SUM;
  int i;
  for(i=0;i<np;i++){
    col[i]=0.;
    nav[i]=0.;
  }
  return 1;
}

int OutputArr::init(double t1, double t2, int npoints, int nc){
  ftype=YDATA;
  strcpy(file,"out");
  ncol=nc;
  np=npoints;
  if(np<2){
    msg_error("OutputArr: too few points (<2)!");
    return 0;
  }
  
  if(c) delete [] c;
  c= new col_desc[ncol];
  if(!c){
    msg_error("OutputArr: MAE\n");
    return 0;
  }
  int i;
  for(i=0;i<ncol;i++){
    c[i].init(np);
    
  }
  ts=t1;
  te=t2;
  dt=(t2-t1)/(npoints-1);
  //col_desc[0].colprop=COLP_STAT; // by default 
  return 1;

}


int  OutputArr::NextPoint(double x, double *y, double *n){
  if(ftype!=XYDATA)return -1;
  double one=1.; 
  pxx.add(&x);
  int i;
  for(i=0;i<ncol;i++){
    pyy[i]->add(y+i);
    if(n)pnn[i]->add(n+i);
    else pnn[i]->add((void *)&one);
  }
  np=pxx.ind;
  return 1;
}


void OutputArr::dealloc(){
  if(ftype==XYDATA){
    int i;
    for(i=0;i<ncol;i++){
      delete pyy[i];
      delete pnn[i];
    }
    delete [] pyy;
    delete [] pnn;
    pxx.clear();
  }
  delete [] c;
  c=NULL;
}


int OutputArr::init(int nc){
  if(c){
    dealloc();
  }

  ftype=XYDATA;
  strcpy(file,"out");
  ncol=nc;
  np=0;

  c= new col_desc[ncol];
  pyy= new PoolP[ncol];
  pnn= new PoolP[ncol];
  if(!c || ! pyy || !pnn){
    msg_error("OutputArr: MAE\n");
    return 0;
  }


  int i;
  for(i=0;i<ncol;i++){
    c[i].colprop=COLP_SUM;
    pyy[i]= new Pool((void **)&c[i].col,(int)sizeof(double));
    pnn[i]= new Pool((void **)&c[i].nav,(int)sizeof(double));
  }

  ts=0;
  te=0;
  dt=0;
  //col_desc[0].colprop=COLP_STAT; // by default 
  return 1;
}



int OutputArr::UpdateColumn(int cn, double F(double)){
  
  if(ftype==XYDATA)return UpdateColumnIR(cn, F);
    

  if(cn<0 || cn>=ncol || c[cn].colprop == COLP_STAT) return -1;
  TableFunction fc(np,c[cn].col), fn(np,c[cn].nav);
  fc.xscale(ts,te);
  fn.xscale(ts,te);
  fc+=F;
  fn+=1;
  return cn;
}
  

int OutputArr::UpdateColumn(int cn, TableFunction *data, TableFunction *xn){

  if(ftype==XYDATA)return UpdateColumnIR(cn,data,xn);

  if(cn<0 || cn>=ncol || c[cn].colprop == COLP_STAT) return -1;
  
  TableFunction fc(np,c[cn].col), fn(np,c[cn].nav);
  fc.xscale(ts,te);
  fn.xscale(ts,te);
  

  fc+=*data;

  if(xn){
    fn+=*xn;
  }
  else{
    fn+=1;
  }

  return cn;
}


int OutputArr::UpdateColumnIR(int cn, double F(double)){
  if(cn<0 || cn>=ncol || c[cn].colprop == COLP_STAT) return -1;
  double *xx=(double *)pxx[0];
  TableFunction fc(np,xx,c[cn].col), fn(np,xx,c[cn].nav);
  fc+=F;
  fn+=1;
  return cn;
}
  

int OutputArr::UpdateColumnIR(int cn, TableFunction *data, TableFunction *xn){
  if(cn<0 || cn>=ncol || c[cn].colprop == COLP_STAT) return -1;
  double *xx=(double *)pxx[0];
  TableFunction fc(np,xx,c[cn].col), fn(np,xx,c[cn].nav);
 
  fc+=*data;
  if(xn){
    fn+=*xn;
  }
  else{
    fn+=1;
  }
  return cn;
}




int OutputArr::PrintToFile(FILE *f){
  double ti,y;
  int i,j,na;
  double *xx=NULL;
  if(ftype==XYDATA)xx=(double *)pxx[0];

  for(i=0;i<np;i++){
    if(ftype==XYDATA)ti=xx[i];
    else ti=ts+i*dt;
    fprintf(f,"%e",ti*arg_coeff);
    for(j=0;j<ncol;j++){
      y=c[j].col[i];
      na=(int)c[j].nav[i];
      if(na)y/=na;
      fprintf(f," %e",y);
    }
    for(j=0;j<ncol;j++){
      na=(int)(c[j].nav[i]+0.00001);
      fprintf(f," %d",na);
    }
    fprintf(f,"\n");
  }
  return 1;
}


class FileCol{
  TableFunction fnc;
public:
  double *arr;
  Pool parr;
  FileCol():parr((void **)&arr,sizeof(double)){
  }
  TableFunction &func(){
    TableFunction f(parr.ind,arr);
    fnc=f;
    return fnc;
  }
};



int OutputArr::Load(){
  FILE *f;

  f=fopen(file,"rt");
  if(!f){
    eprintf("OutputArr:Load: Can't open file '%s'.\n",file);
    return 0;
  }
  char str[1000];
  long pos;
  int lncount=0;

  // skipping header(s)
  int k;
  do{
    pos=ftell(f);
    fgetline(f,str,1000);
    lncount++;
    k=0;
    while(str[k]==' ')k++;
  }while(!isdigit(str[k]) && (str[k]!='.'));

  // counting column number
  char *scan=strtok(str," ");
  int nn=0;
  while(scan){
    scan=strtok(NULL," ");
    nn++;
  }

  // first data column
  int n1=1;
  // first nav column
  int n2=n1+(nn-n1)/2;

  int nact=ncol;
  if(n2-n1<ncol){
    eprintf("OutputArr:Load: warning: Too few columns (%d) in file '%s'.\n"
	    ,n2-n1,file);
    nact=n2-n1;
  }
  else if(n2-n1>ncol){
    eprintf("OutputArr:Load: warning: Too many columns in file '%s'.\n",file);
  }

  int res=1;
  TableFunction *fy, *fn;
  // reading data
  FileCol clx, *cly= new FileCol[nact], *cln= new FileCol[nact];
  if(!cly || !cln){
    msg_error("OutputArr:Load: MAE!\n");
    return 0;
  }
  double ftmp;
  int i;

  fseek(f,pos,SEEK_SET);
  lncount--;

  do{
    fgetline(f,str,1000);
    if(feof(f))break;

    lncount++;
    scan=strtok(str," ");
    // reading x column
    if(!scan ||  !sscanf(scan,"%lf",&ftmp)){
      eprintf("OutputArr:Load: read failure in file"
	      " '%s' (x-column, line #%d)\n",file,lncount);
      res=0;
      break;
    }

    clx.parr.add(&ftmp);
    for(i=0;i<2*nact;i++){
      scan=strtok(NULL," ");
      // reading y (n) column
      if(!scan ||  !sscanf(scan,"%lf",&ftmp)){
	eprintf("OutputArr:Load: read failure in file"
		" '%s' (column #%d of %d, line #%d)\n",file,
		i+2,2*nact,lncount);
	res=0;
	break;
      }
      if(i<nact)cly[i].parr.add(&ftmp);
      else{
	cln[i-nact].parr.add(&ftmp);
	//printf("%f ",ftmp);
      }
    }
    //printf("\n");
  }while(res);

  int n=clx.parr.ind;

  for(i=0;i<nact;i++){
    fy= new TableFunction(n,clx.arr,cly[i].arr);
    fn= new TableFunction(n,clx.arr,cln[i].arr);

    if(!fy || !fn){
      msg_error("OutpuArr:Load: MAE\n");
      return 0;
    }

    (*fy)*=(*fn);
    UpdateColumn(i,fy,fn);

    delete fy;
    delete fn;
  }


  /*
  double *xx=NULL,*yy=NULL, *xn=NULL;

  int n,na,i;
  for(i=0;i<nact;i++){
    fseek(f,pos,SEEK_SET);
    n=ReadTableDirect(f,xx,yy,1,n1+i+1,1);
    fseek(f,pos,SEEK_SET);
    na=ReadTableDirect(f,xn,xn,-1,n2+i+1,1);
    if(n!=na || n<2){
      eprintf("OutputArr:Load: read failure in file '%s'"
	      " (%d, %d)\n",file,n,na);
      res=0;
      break;
    }

    fy= new TableFunction(n,xx,yy);
    fn= new TableFunction(n,xx,xn);

    if(!fy || !fn){
      msg_error("OutpuArr:Load: MAE\n");
      return 0;
    }

    (*fy)*=(*fn);
    UpdateColumn(i,fy,fn);

    delete fy;
    delete fn;
  }
  if(xx) delete [] xx;
  if(yy) delete [] yy;
  if(xn) delete [] xn;
  */

  fclose(f);
  delete [] cly;
  delete [] cln;
  return res;
}

















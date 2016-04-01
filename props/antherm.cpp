# include <stdio.h>
# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "antherm.h"
# include "../interact.h"
# include "../plasmar.h"
# include "statist.h"


int AnTherm::Init(AnalyseShell *parent,int initype){

  if(!Analysator::Init(parent,initype))return 0;


  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&VEL)){
      msg_error("Can't calculate temperature log,\n"
		" trajectory file does not contain Velocity data.\n");
      return 0;
    }
  }


  Set_stop(0);

  if(!Read_param("Electron temperature output step: %lf", &tav))tav=1.;
  if(!Read_param("Ion temperature output step: %lf", &tavi))tavi=30.;
  long aw;
  if(!Read_param("Temperature run-average window: %ld", &aw))aw=10;

  char tmpstr[250];
  if(Read_param("Center-of-mass for components: %s",tmpstr)){
     if(strstr(tmpstr,"separate")){
       printf("Plasma with separated center of masses specified!\n");
       one_center=0;
     }
  }
  else one_center=1;

  if(Read_param("Draw averages: %s",tmpstr)){
     if(strstr(tmpstr,"yes"))draw_av=0;
  }
  else draw_av=1;

  //printf("aw: %ld\n",aw);

  STe.set_window(aw);
  STi.set_window(aw);
  ST.set_window(aw);


  AddIrOutput(2);  // 00 - for Te and T
  AddIrOutput(1);  // 01 for Ti


  t0=psh->rtime-2*tav;
  t0e=psh->rtime;
  t0i=psh->rtime-2*tavi;
  t0iw=psh->rtime;

  eprintf("Init: Temperature log analysator.\n");
  Set_stop(1);
  return 1;
}

int AnTherm::Step(){
  static int cli=0, cle=0;


  if(!(psh->status&PR_IONVEL)){
    if(!cli && psh->rtime-t0i > 2*psh->idt){
      STi.clear();
      cli=1;
    }
  }

  if(!(psh->status&PR_ELCVEL)){
    if(!cle && psh->rtime-t0e > 2*psh->edt){
      STe.clear();
      ST.clear();
      cle=1;
    }
  }



  if(!(psh->status&PR_IONVEL) && !(psh->status&PR_ELCVEL)){
    return 1; // nothing to do
  }

   // setting temperature measurement parameter
  int tmp=psh->gasp->one_center;
  psh->gasp->one_center=one_center;

  psh->gasp->getT();

  if(psh->status&PR_IONVEL){
    cli=0;
    t0i=psh->rtime;
    Ti=psh->gasp->Ti;
    STi.next();
  }

  if(psh->status&PR_ELCVEL){
    cle=0;
    t0e=psh->rtime;
    Te=psh->gasp->Te;
    T=psh->gasp->T;
    ST.next();
    STe.next();
  }

  double y[2];
  //printf("%f %f\n",(psh->rtime-t0),tav);
  if((psh->rtime-t0)>tav){
    //printf("%e %e %e\n",psh->rtime,STe.av(),ST.av());

    //fprintf(f,"%e %e %e\n",psh->rtime,STe.av(),ST.av());
    y[0]=STe.av();
    y[1]=ST.av();
    Out[0].NextPoint(psh->rtime,y);
    t0=psh->rtime;
  }

  if((psh->rtime-t0iw)>tavi){
    //printf("%e %e\n",psh->rtime,STi.av());
    //fprintf(fi,"%e %e\n",psh->rtime,STi.av());
    y[0]=STi.av();
    Out[1].NextPoint(psh->rtime,y);
    t0iw=psh->rtime;
  }

  // restoring temperature measurement parameter
  psh->gasp->one_center=tmp;

  return 1;
}




int  AnTherm::Process(char *name){
  char str[256]="therm.dat";
  strcpy(str,psh->out_dir);
  strcat(str,name);
  strcat(str,"-t.dat");

  
  AssociateOutput(str,0,psh->update);
  UpdateFile(0,"1-time 2-Te 3-T");
  
  //f=Err_fopen(str,"wt");
  //fprintf(f,"#1-time 2-Te 3-T\n");

  strcpy(str,psh->out_dir);
  strcat(str,name);
  strcat(str,"-ti.dat");
  AssociateOutput(str,1,psh->update);
  UpdateFile(1,"1-time 2-Ti");

  //fi=Err_fopen(str,"wt");
  //fprintf(fi,"#1-time 2-Ti\n");



  eprintf("Wrote temperature log...\n");
  //fclose(f);
  //fclose(fi);
  return 1;
}















# include <stdio.h>
# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "antrjlog.h"
# include "../interact.h"
# include "../plasmar.h"


int AnTrjLog::Init(AnalyseShell *parent,int initype){
  
  if(!Analysator::Init(parent,initype))
    return 0;

  eprintf("Init: Animation logger.\n");

  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&COORD)){
      msg_error("Can't produce trajectory log,\n"
		" trajectory file does not contain Coords data.\n");
      return 0;
    }
  }

  //Read_param("Animation file: %s",anim_file);
  Read_paramn(2,"Particle names: %s %s",&name[0],&name[1]);
  Set_stop(0);
  scale=-1;
  Read_param("Animation scale: %lf",&scale);
  Set_stop(1);
  init_file();
  return 1;
}


int AnTrjLog::init_file(){
  char str[1000];
  if(psh->sp_type&PROF_SPLIT)
    sprintf(str,"%s%s%-anim03d.xyz",psh->out_dir,psh->dname,psh->sp_num+1);
  else
    sprintf(str,"%s%s-anim.xyz",psh->out_dir,psh->dname);

  fp=fopen(str,"wt");
  if(!fp){
    msg_error("Animation: Can't open file '%s' for writing!\n",str);
    return 0;
  }
  return 1;
}

int AnTrjLog::Step(){
  if(!(psh->status&PR_IONCOORDS) 
     && !(psh->status&PR_ELCCOORDS))return 1; //nothing to do

  int i;
 
  int ni=psh->gasp->ni; // only for charge-symmetrical case!
  int npart=psh->gasp->n;

  if(!psh->xx_used){
    if(psh->status&PR_IONCOORDS){
      for(i=0;i<ni;i++)psh->gasp->rcell(i);
    }
    if(psh->status&PR_ELCCOORDS){
      for(i=ni;i<npart;i++)psh->gasp->rcell(i);
    }
    //psh->anx=psh->gasp->xx;
    psh->xx_used=1;
  }

  double c= scale >0 ? scale*psh->gasp->L : 1.;

  if(psh->status&PR_ELCCOORDS){
    fprintf(fp,"%d\n",npart);
    fprintf(fp,"t=%g\n", psh->rtime);
    for(i=0;i<npart;i++){
      if(i<ni)
        fprintf(fp,"%s", name[0]);
      else
        fprintf(fp,"%s", name[1]);
      fprintf(fp," %g %g %g\n",c*psh->gasp->xx[i][0],c*psh->gasp->xx[i][1],c*psh->gasp->xx[i][2]);
    }
    fflush(fp);
  }	
  
  
  return 1;
}



int  AnTrjLog::Process(char *name){
  if(!(psh->sp_type&PROF_SPLIT)){
    eprintf("Finishing animation file...\n");
    if(fp)
      fclose(fp);
    fp=NULL;
    eprintf("Done.\n");
  }
  return 1;
}


int AnTrjLog::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    eprintf("Writing animation file for split #%d...\n",nsp);
    if(fp)
      fclose(fp);
    fp=NULL;
    if(!psh->lastsplit)
      init_file();
  }
  return res;
}











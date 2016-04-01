# include <stdio.h>
# include <math.h>

# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "ancurcur.h"
# include "../interact.h"
# include "../plasmar.h"
# include "statist.h"



void AnCurCur::allocate_data(){
  if(vc_ion || vc_tot || vc_cross){
    vc_Ci=new Correlation[3];
    if(!vc_Ci)serror("Current acf: memory allocation error.\n");
  }

  if(vc_elc || vc_tot || vc_cross){
    vc_Ce=new Correlation[3];
    if(!vc_Ce)serror("Current acf: memory allocation error.\n");
    if(vc_tot){
       vc_Ct=new Correlation[3];
       if(!vc_Ct)serror("Current acf: memory allocation error.\n");
    }
    if(vc_cross){
       vc_Cei=new Correlation[3];
       vc_Cie=new Correlation[3];
       if(!vc_Cei || !vc_Cie )serror("Current acf: memory allocation error.\n");
    }
  }
}


// nisteps and nesteps must be known
void AnCurCur::init_data(){
  int i;

  if(vc_ion || vc_tot || vc_cross){
    eprintf("Init: Ion current acf analyser.\n");
    if(vc_type==VC_FFT){
      eprintf("Fourier degree: %d",(int)log2(nisteps)+1);
      vc_ni=1<<((int)log2(nisteps)+1);
      eprintf(", %d points\n",vc_ni);
    }
    else vc_ni=nisteps;

    for(i=0;i<3;i++){
      vc_Ci[i].init(vc_ni,psh->split_length());
      vc_Ci[i].sub_mean=vc_sav;
      vc_Ci[i].insert_begin(nisteps);
    }
  }

  if(vc_elc || vc_tot || vc_cross){
    eprintf("Init: Electron current acf analyser.\n");
    if(vc_type==VC_FFT){
      eprintf("Fourier degree: %d",(int)log2(nesteps)+1);
      vc_ne=1<<((int)log2(nesteps)+1);
      eprintf(", %d points\n",vc_ne);
    }
    else vc_ne=nesteps;

    for(i=0;i<3;i++){
      vc_Ce[i].init(vc_ne,psh->split_length());
      vc_Ce[i].sub_mean=vc_sav;
      vc_Ce[i].insert_begin(nesteps);
    }
    if(vc_tot){
       eprintf("Init: added total current acf analyser.\n");
       for(i=0;i<3;i++){
         vc_Ct[i].init(vc_ne,psh->split_length());
         vc_Ct[i].sub_mean=vc_sav;
         vc_Ct[i].insert_begin(nesteps);
       }
    }
    if(vc_cross){
       eprintf("Init: added cross current acf analyser.\n");
       for(i=0;i<3;i++){
         vc_Cei[i].init(vc_ne,psh->split_length());
         vc_Cei[i].sub_mean=vc_sav;
         vc_Cie[i].init(vc_ne,psh->split_length());
         vc_Cie[i].sub_mean=vc_sav;
       }
    }

  }
}



int AnCurCur::Init(AnalyseShell *parent,int initype){

  if(!Analysator::Init(parent,initype))return 0;


  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&VEL)){
      msg_error("Can't calculate current autocorrelation function,\n"
                "trajectory file does not contain Velocity data.\n");
      return 0;
    }
  }

  Set_stop(1);
  char tmpstr[256];


  Read_param("Ion current acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_ion=1;
  else vc_ion=0;

  Read_param("Electron current acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_elc=1;
  else vc_elc=0;

  Read_param("Total current acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_tot=1;
  else vc_tot=0;

  Read_param("Cross current acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_cross=1;
  else vc_cross=0;

  Set_stop(1);

  Read_param("CACF substract mean: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_sav=1;
  else vc_sav=0;

  Read_param("CACF normalize: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_norm=1;
  else vc_norm=0;


  Read_param("CACF calculation type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))vc_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))vc_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))vc_type=VC_STRAIT;
  else serror("Current acf: unknown calculation type.\n");

  if(!vc_ion && !vc_elc && !vc_tot && !vc_cross){
    eprintf("Warning: nothing to do for current acf!\n");
  }


  allocate_data();
  nisteps=psh->regsteps(PR_IONVEL);
  nesteps=psh->regsteps(PR_ELCVEL);
  init_data();


  if(psh->sp_type&AV_SPLIT){
    eprintf("AnCurCur: warning: split averaging is not implemented!\n");
  }

  return 1;
}


AnCurCur::~AnCurCur(){
  if(vc_ion || vc_tot)delete[] vc_Ci;
  if(vc_elc || vc_tot)delete[] vc_Ce;
  if(vc_tot)delete[] vc_Ct;
  if(vc_cross){
     delete[] vc_Cei;
     delete[] vc_Cie;
  }
}


int AnCurCur::Step(){
  if(!(psh->status&PR_IONVEL)
     && !(psh->status&PR_ELCVEL))return 1; //nothing to do

  int nion=psh->gasp->ni;
  int nelec=psh->gasp->ne;

  Vector_3 sumvi(0,0,0), sumve(0,0,0);
  int i,j;
  if((psh->status&PR_IONVEL) && (vc_ion||vc_tot)){
   for(i=0;i<nion;i++){
     sumvi+=psh->anv[i];
   }
   for(j=0;j<3;j++){
     vc_Ci[j].insert_next(sumvi[j]);
   }
  }

  if((psh->status&PR_ELCVEL) && (vc_elc||vc_tot)){
   for(i=0;i<nelec;i++){
     sumve+=psh->anv[nion+i];
   }
   for(j=0;j<3;j++){
     vc_Ce[j].insert_next(sumve[j]);
     if(vc_tot)vc_Ct[j].insert_next(sumve[j]); // then we will call prepare_total to make current
   }
  }
  return 1;
}




int AnCurCur::output_vcorr(char *fileform,Correlation* Corr){
  int i,j;

  char filename[250];
  double *acf_data;

  for(i=1;i<3;i++){
    for(j=0;j<Corr[0].n;j++){
       Corr[0].arr[j]+=Corr[i].arr[j];
    }
  }
  if(vc_norm)Corr[0].normalize();
  acf_data=Corr[0].arr;

  if(vc_type==VC_DIRECT ){
    sprintf(filename,fileform,"d");
  }
  else if(vc_type==VC_STRAIT){
    sprintf(filename,fileform,"s");
  }
  else{
    sprintf(filename,fileform,"f");
  }
  

  double tmax=Corr[0].maxtime;

  TableFunction acf(Corr[0].n,acf_data);
  acf.xscale(0.,tmax);

  TableFunction cur(Corr[0].n,Corr[0].ini);
  cur.xscale(0.,tmax);

  if(vc_type==VC_FFT)tmax/=2;

  double dt=fmax(0.1,tmax/Corr[0].n);

  int npo=(int)(tmax/dt)+1;
  if(npo<2)npo=2;
  //int npo=Corr[0].n;

  if(!ExistsOutput(0))AddOutput(0,tmax,npo,3);
  else Out[0].init(0,tmax,npo,1);

  AssociateOutput(filename,0,psh->update);
  Out[0].UpdateColumn(0,&acf);
  acf.integrate(1./(4*M_PI));
  Out[0].UpdateColumn(1,&acf);

  
  Out[0].UpdateColumn(2,&cur);

  UpdateFile(0,"1-t*Wpe  2-<j(t)j(0)>/<j(0)j(0)>  3-sigma(t) 4-current");
  Out[0].dealloc();

  return 1;
}


// making total current
void AnCurCur::prepare_total(){
  TableFunction tmpf;
  double q=psh->gasp->q;
  int i;
  for(i=0;i<3;i++){
    tmpf<<vc_Ci[i].fini;
    tmpf*=-q;
    vc_Ct[i].fini+=tmpf;
    tmpf<<vc_Ci[i].finii;
    tmpf*=-q;
    vc_Ct[i].finii+=tmpf;
  }
}

void AnCurCur::insert_end(){
  int j;
  for(j=0;j<3;j++){
    if(vc_elc || vc_tot || vc_cross)vc_Ce[j].insert_end();
    if(vc_ion || vc_tot || vc_cross)vc_Ci[j].insert_end();
  }
}


void AnCurCur::do_correlations(){
  int j;

  if(vc_type==VC_DIRECT){
    eprintf("Direct correlation ...\n");
    for(j=0;j<3;j++){
      //Corr[j].insert_end();
      if(vc_cross){
        CrossCorrDirect(&vc_Cei[j],&vc_Cie[j],&vc_Ce[j],&vc_Ci[j]);
      }
      if(vc_elc || vc_tot)vc_Ce[j].direct();
      if(vc_ion || vc_tot)vc_Ci[j].direct();
      if(vc_tot)vc_Ct[j].direct();
    }
  }
  else if(vc_type==VC_STRAIT){
    eprintf("Strait correlation ...\n");
    for(j=0;j<3;j++){
      //Corr[j].insert_end();
      if(vc_cross){
        CrossCorrStrait(&vc_Cei[j],&vc_Cie[j],&vc_Ce[j],&vc_Ci[j]);
      }
      if(vc_elc || vc_tot)vc_Ce[j].strait();
      if(vc_ion || vc_tot)vc_Ci[j].strait();
      if(vc_tot)vc_Ct[j].strait();
    }
  }
  else if(vc_type==VC_FFT){
    eprintf("FFT correlation ...\n");
    for(j=0;j<3;j++){
      //Corr[j].insert_end();
      if(vc_cross){
        CrossCorrFFT(&vc_Cei[j],&vc_Cie[j],&vc_Ce[j],&vc_Ci[j]);
      }
      if(vc_elc || vc_tot)vc_Ce[j].calculate();
      if(vc_ion || vc_tot)vc_Ci[j].calculate();
      if(vc_tot)vc_Ct[j].calculate();
    }
  }
}



int  AnCurCur::Process(char *name){
  int res=1;
  if(!(psh->sp_type&PROF_SPLIT)){
    char str[256]="vcorr.dat";

    insert_end();
    if(vc_tot)prepare_total();
    do_correlations();

    if(vc_tot){
      // calculating total current
      sprintf(str,"%s%s-jj%%s.dat",psh->out_dir,name);
      eprintf("Writing total current correlation...\n");
      res&=output_vcorr(str,vc_Ct);
    }
    if(vc_ion){
      sprintf(str,"%s%s-jji%%s.dat",psh->out_dir,name);
      eprintf("Writing ion current correlation...\n");
      res&=output_vcorr(str,vc_Ci);
    }
    if(vc_elc){
      sprintf(str,"%s%s-jje%%s.dat",psh->out_dir,name);
      eprintf("Writing electron current correlation...\n");
      res&=output_vcorr(str,vc_Ce);
    }
    if(vc_cross){
      sprintf(str,"%s%s-jjei%%s.dat",psh->out_dir,name);
      eprintf("Writing e-i current correlation...\n");
      res&=output_vcorr(str,vc_Cei);
      sprintf(str,"%s%s-jjie%%s.dat",psh->out_dir,name);
      eprintf("Writing i-e current correlation...\n");
      res&=output_vcorr(str,vc_Cie);
    }

    eprintf("Done.\n");
  }
  return res;
}


int  AnCurCur::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    char str[256]="vcorr.dat";

    insert_end();
    if(vc_tot)prepare_total();
    do_correlations();

    if(vc_tot){
      // calculating total current
      sprintf(str,"%s%s-jj%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing total current correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Ct);
    }
    if(vc_ion){
      sprintf(str,"%s%s-jji%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing total current correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Ci);
    }
    if(vc_elc){
      sprintf(str,"%s%s-jje%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing total current correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Ce);
    }
    if(vc_cross){
      sprintf(str,"%s%s-jjei%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing total current correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Cei);
      sprintf(str,"%s%s-jj%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing total current correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Cie);
    }
    if(!psh->lastsplit){
       eprintf("Reinitializing for new split...\n");
       nisteps=psh->regsteps_next(PR_IONVEL);
       nesteps=psh->regsteps_next(PR_ELCVEL);
       init_data();
    }
    eprintf("Done.\n");
  }
  return res;
}













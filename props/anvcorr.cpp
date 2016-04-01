# include <stdio.h>
# include <math.h>

# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "anvcorr.h"
# include "../interact.h"
# include "../plasmar.h"
# include "statist.h"


int AnVcorr::Init(AnalyseShell *parent,int initype){
  
  if(!Analysator::Init(parent,initype))return 0;


  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&VEL)){
      msg_error("Can't calculate velocity distribution function,\n"
			 "trajectory file does not contain Velocity data.\n");
      return 0;
    }
  }

  Set_stop(0);
  int i;
  char tmpstr[256];

  Read_param("Ion velocity acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_ion=1;
  else vc_ion=0;

  Read_param("Electron velocity acf: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_elc=1;
  else vc_elc=0;

  Set_stop(1);

  Read_param("Substract mean: %s",tmpstr);
  if(strstr(tmpstr,"y"))vc_sav=1;
  else vc_sav=0;


  Read_param("Calculation type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))vc_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))vc_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))vc_type=VC_STRAIT;
  else serror("Velocity acf: unknown calculation type: %s.\n",tmpstr);

  if(!vc_ion && !vc_elc){
    eprintf("Warning: nothing to do for velocity acf!\n");
  }

  int nion, nelec;
  
  if(vc_ion){
    nisteps=psh->regsteps(PR_IONVEL);
    nion=psh->gasp->ni;
    eprintf("Init: Ion velocity acf analyser.\n");
    if(vc_type==VC_FFT){
      eprintf("Fourier degree: %d",(int)log2(nisteps)+1);
      vc_ni=1<<((int)log2(nisteps)+1);
      eprintf(", %d points\n",vc_ni);
    }
    else vc_ni=nisteps;

    vc_Ci=new Correlation[3*nion];
    if(!vc_Ci)serror("Velocity acf: memory allocation error.\n");
    for(i=0;i<3*nion;i++){
      vc_Ci[i].init(vc_ni,psh->split_length());
      vc_Ci[i].sub_mean=vc_sav;
      vc_Ci[i].insert_begin(nisteps);
    }
  }

  if(vc_elc){
    nesteps=psh->regsteps(PR_ELCVEL);
    nelec=psh->gasp->ne;
    eprintf("Init: Electron velocity acf analyser.\n");
    if(vc_type==VC_FFT){
      eprintf("Fourier degree: %d",(int)log2(nesteps)+1);
      vc_ne=1<<((int)log2(nesteps)+1);
      eprintf(", %d points\n",vc_ne);
    }
    else vc_ne=nesteps;
    vc_Ce=new Correlation[3*nelec];
    if(!vc_Ce)serror("Velocity acf: memory allocation error.\n");
    for(i=0;i<3*nelec;i++){
      vc_Ce[i].init(vc_ne,psh->split_length());
      vc_Ce[i].sub_mean=vc_sav;
      vc_Ce[i].insert_begin(nesteps);
    }
  } 

  if(psh->sp_type&AV_SPLIT){
    eprintf("AnVcorr: warning: split averaging is not implemented!\n");
  }
  
  return 1;
}


AnVcorr::~AnVcorr(){
  if(vc_ion)delete[] vc_Ci;
  if(vc_elc)delete[] vc_Ce;
}


int AnVcorr::Step(){
  if(!(psh->status&PR_IONVEL) 
     && !(psh->status&PR_ELCVEL))return 1; //nothing to do

  int nion=psh->gasp->ni;
  int nelec=psh->gasp->ne;

  int i,j;
  if((psh->status&PR_IONVEL) && vc_ion){
   for(i=0;i<nion;i++){
    for(j=0;j<3;j++){
     vc_Ci[3*i+j].insert_next(psh->anv[i][j]);
    }
   }
  }
  if((psh->status&PR_ELCVEL) && vc_elc){
   for(i=0;i<nelec;i++){
    for(j=0;j<3;j++){
     vc_Ce[3*i+j].insert_next(psh->anv[nion+i][j]);
    }
   }
  }
  return 1;
}




int AnVcorr::output_vcorr(char *fileform,Correlation* Corr,int npart){
  int i,j;

  char filename[250];
  double *acf_data;
 
  if(vc_type==VC_DIRECT){
   sprintf(filename,fileform,"d");
   eprintf("Direct correlation ...\n");
   for(i=0;i<npart;i++){
    //printf("Direct correlation #%d...\n",i+1);
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
     Corr[3*i+j].direct();
    }
   }
   for(j=0;j<Corr[0].n;j++){
    //Corr[0].ini[j]=Corr[0].arr[j]*Corr[0].arr[j];
    for(i=1;i<3*npart;i++){
     Corr[0].arr[j]+=Corr[i].arr[j];
     //Corr[0].ini[j]+=Corr[i].arr[j]*Corr[i].arr[j]; // average square
    }
    Corr[0].arr[j]/=(double)(3*npart);
    //Corr[0].ini[j]/=(double)(3*npart);
    //Corr[0].ini[j]=Corr[0].ini[j]-Corr[0].arr[j]*Corr[0].arr[j];
   }
   Corr[0].normalize();
   //x*=x;
   //for(j=0;j<Corr[0].n;j++)Corr[0].ini[j]/=x;
   acf_data=Corr[0].arr;
  }
  else if(vc_type==VC_STRAIT){
   sprintf(filename,fileform,"s");
   printf("Strait correlation ...\n");
   for(i=0;i<npart;i++){
    //printf("Strait correlation #%d...\n",i+1);
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
     Corr[3*i+j].strait();
    }
   }
   for(j=0;j<Corr[0].n;j++){
    //Corr[0].ini[j]=Corr[0].arr[j]*Corr[0].arr[j];
    for(i=1;i<3*npart;i++){
     Corr[0].arr[j]+=Corr[i].arr[j];
     //Corr[0].ini[j]+=Corr[i].arr[j]*Corr[i].arr[j]; // average square
    }
    Corr[0].arr[j]/=(double)(3*npart);
    //Corr[0].ini[j]/=(double)(3*npart);
    //Corr[0].ini[j]=Corr[0].ini[j]-Corr[0].arr[j]*Corr[0].arr[j];
   }
   Corr[0].normalize();
   //x*=x;
   //for(j=0;j<Corr[0].n;j++)Corr[0].ini[j]/=x;
   acf_data=Corr[0].arr;
  }
  else{
   sprintf(filename,fileform,"f");

   for(i=0;i<npart;i++){
    for(j=0;j<3;j++){
     Corr[3*i+j].insert_end();
    }
   }

   printf("Velocity auto corelations: performing  FFT...\n");
   acf_data=Calculate_forth(Corr[0].n,Corr,3*npart);
   printf("Inverse FFT...\n");
   Calculate_back(acf_data,Corr[0].n);
   double c=acf_data[0];
   for(i=0;i<Corr[0].n;i++)acf_data[i]=acf_data[2*i]/c;
  }
  TableFunction acf(Corr[0].n,acf_data);

  double tmax=Corr[0].maxtime;
  acf.xscale(0.,tmax);

  if(vc_type==VC_FFT)tmax/=2;

  double dt=fmax(0.1,tmax/Corr[0].n);

  int npo=(int)(tmax/dt)+1;
  if(npo<2)npo=2;
  

  if(!ExistsOutput(0))AddOutput(0,tmax,npo,1);
  else Out[0].init(0,tmax,npo,1);

  AssociateOutput(filename,0,psh->update);
  Out[0].UpdateColumn(0,&acf);
  UpdateFile(0,"1-t*Wpe  2-<v(t)v(0)>/<v(0)v(0)>");
  Out[0].dealloc();


# if 0
  //Corr[0].fini.xscale(0.,tmax);
  FILE *f1=Err_fopen(filename,"wt");
  fprintf(f1,"#1-t*Wpe  2-<v(t)v(0)>/<v(0)v(0)>\n");
  double dt=fmax(0.1,tmax/Corr[0].n);
  double t;
  for(t=0.;t<tmax /*&& t< 150.*/;t+=dt){
   
    fprintf(f1,"%f  %f",t,acf(t));
    //if(vc_type==VC_STRAIT || vc_type==VC_DIRECT){
    // fprintf(f1," %f",Corr[0].fini(t));
    //}
    fprintf(f1,"\n");
  }
  fclose(f1);

# endif

  return 1;
}

  

int  AnVcorr::Process(char *name){
  int res=1;
  if(!(psh->sp_type&PROF_SPLIT)){
    char str[256]="vcorr.dat";

    if(vc_ion){
      sprintf(str,"%s%s-vvi%%s.dat",psh->out_dir,name);
      eprintf("Writing ion velocity correlation...\n");
      res&=output_vcorr(str,vc_Ci,psh->gasp->ni);
     
    }
    if(vc_elc){
      sprintf(str,"%s%s-vve%%s.dat",psh->out_dir,name);
      eprintf("Writing electron velocity correlation...\n");
      res&=output_vcorr(str,vc_Ce,psh->gasp->ne);
    }
    eprintf("Done.\n");
  }
  //if(vc_ion)delete[] vc_Ci;
  //if(vc_elc)delete[] vc_Ce;
  
  return res;
}


int  AnVcorr::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    int i;
    char str[256]="vcorr.dat";

    if(vc_ion){
      sprintf(str,"%s%s-vvi%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing ion velocity correlation for split #%d...\n",nsp);
     

      res&=output_vcorr(str,vc_Ci,psh->gasp->ni);

      if(!psh->lastsplit){
	nisteps=psh->regsteps_next(PR_IONVEL);
	eprintf("ReInit: Ion velocity acf analyser %d.\n",nisteps);
	if(vc_type==VC_FFT){
	  eprintf("Fourier degree: %d",(int)log2(nisteps)+1);
	  vc_ni=1<<((int)log2(nisteps)+1);
	  eprintf(", %d points\n",vc_ni);
	}
	else vc_ni=nisteps;
	
	for(i=0;i<3*psh->gasp->ni;i++){
	  vc_Ci[i].init(vc_ni,psh->split_length_next());
	  vc_Ci[i].sub_mean=vc_sav;
	  vc_Ci[i].insert_begin(nisteps);
	}

      }
    } 
    if(vc_elc){
      sprintf(str,"%s%s-vve%%s.%03d",psh->out_dir,psh->dname,nsp);
      eprintf("Writing electron velocity correlation for split #%d...\n",nsp);
      res&=output_vcorr(str,vc_Ce,psh->gasp->ne);
     
      if(!psh->lastsplit){
	nesteps=psh->regsteps_next(PR_ELCVEL);
	eprintf("ReInit: Electron velocity acf analyser %d.\n",nesteps);
	if(vc_type==VC_FFT){
	  eprintf("Fourier degree: %d",(int)log2(nesteps)+1);
	  vc_ne=1<<((int)log2(nesteps)+1);
	  eprintf(", %d points\n",vc_ne);
	}
	else vc_ne=nesteps;
	
	for(i=0;i<3*psh->gasp->ne;i++){
	  vc_Ce[i].init(vc_ne,psh->split_length_next());
	  vc_Ce[i].sub_mean=vc_sav;
	  vc_Ce[i].insert_begin(nesteps);
	}
      }
    }
    eprintf("Done.\n");
  }
  return res;
}













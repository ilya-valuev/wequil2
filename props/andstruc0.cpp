# include <stdio.h>
# include <math.h>

# include "common.h"
# include "rparam.h"

# ifdef UNIX
# include "../plasma.h"
# include "andstruc.h"
# include "../interact.h"
# include "../plasmar.h"
# include "../statist.h"
# else
# include "plasma.h"
# include "andstruc.h"
# include "interact.h"
# include "plasmar.h"
# include "statist.h"
# endif

# include "four.h"

# define MAXLEV 5
# define MAXSPLITS 200


char *AnDsymb[]={
  "Sii","See","Sie","Sei","S","Z"
};

unsigned arrp=0;

int _splitsq(unsigned n, int lev, unsigned *arr,unsigned limit){
  float lsq=n;
  lsq=sqrt(lsq);  
  unsigned rr,rest;
  rr=(unsigned)lsq;

  rest=n-rr*rr;
  if(lev==0 && rest !=0)return 0;

  int nsp=0,i;
  do{
    arr[arrp*MAXLEV+lev]=rr;
    if(rest==0){
      for(i=lev-1;i>=0;i--)arr[arrp*MAXLEV+i]=0;
      if(arrp>=limit-1)printf("splitsq: Array limit reached!\n");
      else arrp++;
      for(i=MAXLEV-1;i>=lev;i--)arr[arrp*MAXLEV+i]=arr[(arrp-1)*MAXLEV+i];
      nsp++;
      if(lev==0 || rr==0)return nsp;
    }
    else nsp+=_splitsq(rest,lev-1,arr,limit);

    if(rr==0)break;
    rr--;
    rest=n-rr*rr;

  }while(1/*rr!=0 && rest<(rr+1)*(rr+1)*/);

  return nsp;
} 



int splitsq(unsigned n, int lev,unsigned *arr,unsigned limit=MAXSPLITS){
  arrp=0;
  if(lev>MAXLEV){
    printf("spltsq: Level is reduced to MAXLEV=%d\n",MAXLEV);
    lev=MAXLEV;
  }
  return _splitsq(n,lev-1,arr,limit);
}    


double random2(void){
 return 1.-2.*((double)rand())/RAND_MAX;
}



int AnDstruc::Init(AnalyseShell *parent,int initype){
  
  if(!Analysator::Init(parent,initype))return 0;

  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&COORD)){
      msg_error("Can't calculate dynamical structure factor,\n"
			 "trajectory file does not contain Coords data.\n");
      return 0;
    }
  }

  eprintf("Init: Dynamical structure factor analyser.\n");

  unsigned int arr[MAXLEV*MAXSPLITS];
  Set_stop(1);

  char tmpstr[200];
  Read_param("DSF calculation type: %s",tmpstr);
  if(strstr(tmpstr,"fft"))nd_type=VC_FFT;
  else if(strstr(tmpstr,"direct"))nd_type=VC_DIRECT;
  else if(strstr(tmpstr,"strait"))nd_type=VC_STRAIT;
  else serror("DSF: unknown calculation type.\n");

  wr_tc=1;
  Read_param("DSF write time correlations: %s",tmpstr);
  if(strstr(tmpstr,"n"))wr_tc=0;

  wr_fc=0;
  Read_param("DSF write frequency correlations: %s",tmpstr);
  if(strstr(tmpstr,"y")){
    wr_fc=1;
    Read_param("DSF frequency interval, step:>",tmpstr);
    if(strstr(tmpstr,"correct"))f_corr_int=1;
    else{
      if(strstr(tmpstr,"default")){
	f_corr_int=-1;
	if(sscanf(tmpstr,"%lf, %lf",&f_ws, &f_we)!=2){
	  msg_error("Can't read DSF frequency interval (default step)!\n");
	  return 0;
	}
      }
      else{
	f_corr_int=0;
	if(sscanf(tmpstr,"%lf, %lf, %lf",&f_ws, &f_we, &f_dw)!=3){
	  msg_error("Can't read DSF frequency interval and step!\n");
	  return 0;
	}
      }
      if(f_ws>=f_we){
	msg_error("Incorrect DSF frequency interval (%f, %f)!\n",f_ws,f_we);
	return 0;
      }
    }
    Read_param("DSF integration time: %s",tmpstr);
    if(strstr(tmpstr,"half_split"))f_itime=-1.;
    else{
      if(!sscanf(tmpstr,"%lf",&f_itime)){
	msg_error("Can't read DSF integration time!\n");
	return 0;
      }
    }
    Read_param("DSF frequency output:>",tmpstr);
    char str[50];
    f_outn=0;
    int i;
    for(i=0;i<6;i++){
      f_out[i]=FOUT_NONE;
      sprintf(str,"%s%s",AnDsymb[i],"r");
      if(strstr(tmpstr,str)){
	f_out[i]=FOUT_R;
	f_outn++;
      }
      sprintf(str,"%s%s",AnDsymb[i],"i");
      if(strstr(tmpstr,str)){
	f_out[i]=FOUT_I;
	f_outn++;
      }
      sprintf(str,"%s%s",AnDsymb[i],"m");
      if(strstr(tmpstr,str)){
	f_out[i]=FOUT_M;
	f_outn++;
      }
    }
    if(f_outn==0){
      eprintf("Warning: no DSF frequency output specified,\n"
	      "switching frequency DSF off...\n");
      wr_fc=0;
    }
    Set_stop(0);
    if(!Read_param("DSF swoothing window (Wp): %lf",&out_w))out_w=-1;
    Set_stop(1);
  }



  nisteps=psh->regsteps(PR_IONCOORDS);
  nesteps=psh->regsteps(PR_ELCCOORDS);
  nsteps=nesteps;
  //nisteps=nesteps;


  if(nd_type==VC_FFT){
    eprintf("DSF: Electron Fourier degree: %d",(int)log2(nsteps)+1);
    nd_n=1<<((int)log2(nsteps)+1);
    eprintf(", %d points\n",nd_n);

    eprintf("DSF: Ion Fourier degree: %d",(int)log2(nisteps)+1);
    nd_ni=1<<((int)log2(nisteps)+1);
    eprintf(", %d points\n",nd_ni);


  }
  else{
    nd_n=nsteps;
    nd_ni=nisteps;
  }
  nd_ni=nd_n; // assume for simplicity, linear interpolation

  // distribution renormalization factors (must be equal!!!)
  //factor_i=sqrt((float)nd_ni/nisteps);
  factor_e=1; //1./psh->split_length();  //1; //sqrt((float)nd_n/nesteps);
  factor_i=factor_e;

  Set_stop(0);

  int i,l;
  float kval;
  if(Read_param("DSF Free K: %lf %d", &kval, &nd_nvect)==2){
    eprintf("DSF: Warning: using not allowed k-directions! \n");
    nd_kvect= new Vector_3[nd_nvect];
    if(!nd_kvect)serror("DSF: MAE\n");
    for(i=0;i<nd_nvect;i++){
      nd_kvect[i]=Vector_3(random2(),random2(),random2());
      nd_kvect[i].normalize();
      nd_kvect[i]*=kval*2*M_PI/psh->gasp->L;
    }
  }
  else{

    Set_stop(1);
    int ksq;
    Read_param("DSF K square: %d",&ksq);
    //Read_param("Frequency interval, steps: %f, %f, %d",&ds_f0,&ds_f1,&ds_nf);

    int nk;
    nk=splitsq(ksq,3,arr);
    eprintf("%d components\n",nk);



    nd_kvect= new Vector_3[8*nk];
    if(!nd_kvect)serror("DSF: MAE\n");

    nd_nvect=0;
    Vector_3 k;

    for(l=0;l<nk;l++){ //all possible vectors
      for(i=0;i<3;i++){
	k[i]=(float)arr[MAXLEV*l+i];
      }
      k*=2*M_PI/psh->gasp->L;

      // determination of nvect
      int e;
      int m=0;
      while(k[m]==0. && m<3)m++;

      for(e=0;e<4;e++){ // averaging through 8 possible directions

       nd_kvect[nd_nvect++]=k;
       if(e!=0)k[m]=-k[m];

       m=(m+1)%3;
       if(k[m]==0.){
         e+=2;
         m=(m+1)%3;
         if(k[m]==0.)break;
       }
       k[m]=-k[m];
      }
    }

  }

  eprintf("%d vectors\n",nd_nvect);
  if(nd_nvect==0)serror("!!!!!!???????\n");

  int m;
  for(m=0;m<NCORR;m++){
   nd_r[m] = new Correlation[nd_nvect];
   if(!nd_r[m])serror("DSF: MAE\n");


   for(i=0;i<nd_nvect;i++){
    if(m==0){
      nd_r[m][i].init(nd_ni,psh->split_length(),CORR_COMPLEX);
      nd_r[m][i].insert_begin(nisteps,factor_i);
    }
    else{
      nd_r[m][i].init(nd_n,psh->split_length(),CORR_COMPLEX);
      // only for ee and ii
      if(m<2 || m>3)nd_r[m][i].insert_begin(nsteps,factor_e);
    }
   }
  }

  cur_isr=cur_isi=0.;

  if(psh->sp_type&AV_SPLIT){
    eprintf("AnDstruc: warning: split averaging is not implemented!\n");
  }

  return 1;

}



AnDstruc::~AnDstruc(){
  if(nd_nvect){
    int m;
    for(m=0;m<NCORR;m++){
      delete [] nd_r[m];
    }
  }
}


int AnDstruc::Step(){
  if(!(psh->status&PR_IONCOORDS) 
     && !(psh->status&PR_ELCCOORDS))return 1; //nothing to do

  int nion=psh->gasp->ni;
  int nelec=psh->gasp->ne;
  int npart=nion+nelec;

  int i;
  if(!psh->xx_used){
    if(psh->status&PR_IONCOORDS){
      for(i=0;i<nion;i++)psh->gasp->rcell(i);
    }
    if(psh->status&PR_ELCCOORDS){
      for(i=nion;i<npart;i++)psh->gasp->rcell(i);
    }
    //psh->gasp->anx=psh->gasp->xx;
    psh->xx_used=1;
  }


  int j,l;
  Vector_3 r,k;
  float sr=0,si=0;

  for(l=0;l<nd_nvect;l++){ // assembling data for  all possible vectors
    k=nd_kvect[l];

    if(psh->status&PR_IONCOORDS){
      sr=si=0.;
      for(j=0;j<nion;j++){ // averaging through nn2 particles
	r=psh->anx[j];
	
	sr+=cos(r*k);
	si+=sin(r*k);
      }
      nd_r[0][l].insert_next(sr,si);
      cur_isr=sr;
      cur_isi=si;
    }
    //nd_r[0][l].insert_next(cur_isr,cur_isi);

    if(psh->status&PR_ELCCOORDS){
      sr=si=0.;
      for(j=nion;j<npart;j++){ // averaging through nn2 particles
	r=psh->anx[j];

	sr+=cos(r*k);
	si+=sin(r*k);
      }
    }

    nd_r[1][l].insert_next(sr,si);
    nd_r[4][l].insert_next(cur_isr-sr,cur_isi-si);
  }

  return 1;
}



float make_mod(float x,float y){
  return sqrt(x*x + y*y);
}



int AnDstruc::fileout(char *fileform){

  int i;

  for(i=0;i<nd_nvect;i++){
   nd_r[0][i].insert_end();
   nd_r[1][i].insert_end();
  }

  char filename[200];
  char tch;


  if(nd_type==VC_DIRECT){
   tch='d';


   eprintf("DSF: Direct correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     eprintf("Vector #%d\n",i+1);
     CrossCorrDirect(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].direct();
   }
  }
  else if(nd_type==VC_STRAIT){
   tch='s';

   eprintf("DSF: Strait correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     eprintf("Vector #%d\n",i+1);
     CrossCorrStrait(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].strait();
   }
  }
  else{ //FFT
   tch='f';

   eprintf("DSF: FFT correlation ...\n");
   for(i=0;i<nd_nvect;i++){
     eprintf("Vector #%d\n",i+1);
     CrossCorrFFT(&nd_r[2][i],&nd_r[3][i],&nd_r[0][i],&nd_r[1][i]);
     nd_r[4][i].calculate();
   }
  }

  // averaging over vectors

  int j,m;

  for(m=0;m<NCORR;m++){
   int n=nd_r[m][0].n;

   for(j=0;j<n;j++){


    for(i=1;i<nd_nvect;i++){
     nd_r[m][0].arr[j]+=nd_r[m][i].arr[j]; // real part
     nd_r[m][0].arri[j]+=nd_r[m][i].arri[j]; // imaginary part
     if(nd_type==VC_FFT){
      nd_r[m][0].arr[n+j]+=nd_r[m][i].arr[n+j]; // fft^2
      // only for cross-correlations
      if(m>=2)nd_r[m][0].arri[n+j]+=nd_r[m][i].arri[n+j]; // imaginary of fft1 mal fft2*
     }
    }
    nd_r[m][0].arr[j]/=nd_nvect;  //*n;
    nd_r[m][0].arri[j]/=nd_nvect;

    if(nd_type==VC_FFT){
     nd_r[m][0].arr[n+j]/=nd_nvect;
     if(m>=2)nd_r[m][0].arri[n+j]/=nd_nvect;
    }
   }
  }

  float tmax=nd_r[1][0].maxtime; // tmax and n from <ee>
  float dt=fmax(0.02,tmax/nd_r[1][0].n);

  int npo=(int)(tmax/dt)+1;
  if(npo<2)npo=2;

  if(!ExistsOutput(0))AddOutput(0,tmax,npo,12);
  else Out[0].init(0,tmax,npo,12);

  char tst[5]="dt";
  if(wr_tc){
    tst[0]=tch;
    sprintf(filename,fileform,tst);
    eprintf("Updating time correlations file...\n");
    AssociateOutput(filename,0,psh->update);
  }

  TableFunction Sr, Si; // full correlations S=<ii>+<ee>-<ei>-<ie>
  Sr<< nd_r[0][0].func_r();
  Sr+= nd_r[1][0].func_r();
  Sr-= nd_r[2][0].func_r();
  Sr-= nd_r[3][0].func_r();

  Si<< nd_r[0][0].func_i();
  Si+= nd_r[1][0].func_i();
  Si-= nd_r[2][0].func_i();
  Si-= nd_r[3][0].func_i();


  Out[0].UpdateColumn(0,&nd_r[0][0].func_r());  // Re(Sii)
  Out[0].UpdateColumn(1,&nd_r[0][0].func_i());  // Im(Sii)
  Out[0].UpdateColumn(2,&nd_r[1][0].func_r());  // Re(See)
  Out[0].UpdateColumn(3,&nd_r[1][0].func_i());  // Im(See)
  Out[0].UpdateColumn(4,&nd_r[2][0].func_r());  // Re(Sei)
  Out[0].UpdateColumn(5,&nd_r[2][0].func_i());  // Im(Sei)
  Out[0].UpdateColumn(6,&nd_r[3][0].func_r());  // Re(Sie)
  Out[0].UpdateColumn(7,&nd_r[3][0].func_i());  // Im(Sie)

  Out[0].UpdateColumn(8,&Sr);  // Re(S)
  Out[0].UpdateColumn(9,&Si);  // Im(S)

  Out[0].UpdateColumn(10,&nd_r[4][0].func_r());  // Re(S), direct sum
  Out[0].UpdateColumn(11,&nd_r[4][0].func_i());  // Im(S)


  if(wr_tc){
    UpdateFile(0,"1-t 2-Re(Sii) 3-Im(Sii) 4-Re(See) 5-Im(See) 6-Re(Sie) "
	         "7-Im(Sie) 8-Re(Sei) 9-Im(Sei) 10-Re(S) 11-Im(S)");
  }
  if(wr_fc){
    // Fourrier
    // tmax and n from <ee> or f_itime
    if(f_itime<0)tmax=nd_r[1][0].maxtime/2;
    else tmax=fmin(nd_r[1][0].maxtime,f_itime);

    dt=fmax(0.02,nd_r[1][0].maxtime/nd_r[1][0].n);

    int npoc;
    npoc=(int)(tmax/dt)+1;
    if(npoc<2)npoc=2;


    int npof=npoc;
    if(f_corr_int==1){  // adjusting Fourrier interval and step
      f_ws=0;
      f_dw=2*M_PI/tmax;
      f_we=f_ws+f_dw*(npoc-1);
    }
    else if(f_corr_int==0){
      npof=(int) ((f_we-f_ws)/f_dw) +1;
    }
    else{
      f_dw=(f_we-f_ws)/(npoc-1);
    }

    if(!ExistsOutput(1))AddOutput(f_ws,f_we,npof,f_outn);
    else Out[1].init(f_ws,f_we,npof,f_outn);

    tst[0]=tch;
    tst[1]='f';
    sprintf(filename,fileform,tst);
    eprintf("Updating DSF file...\n");
    if(!wr_tc)AssociateOutput(filename,1,psh->update);
    else AssociateOutput(filename,1,0);


    int dim=fmax(npoc,npof);
    float *arr= new float[2*dim], *tmp;
    if(!arr){
      msg_error("AnDstruc:fileout: MAE!\n");
      return 0;
    }

    //TableFunction Fr(dim,arr), Fi(dim,arr+dim);
    TableFunction Fr(npof,arr), Fi(npof,arr+dim);
    Fr.xscale(f_ws,f_we);
    Fi.xscale(f_ws,f_we);
    // determining smoothing window
    int iwin=1;
    if(out_w>0){
      iwin=out_w/f_dw+1;
    }

    eprintf("Frequency transform...\n");
    char str[500]="1-w", str2[250];
    int k=0;

    for(m=0;m<6;m++){
      if(f_out[m]==FOUT_NONE)continue;

      eprintf("Columns %d and %d...\n",2*m+2, 2*m+3);
      for(i=0;i<npoc;i++){
	arr[2*i  ]=Out[0].c[2*m  ].col[i];
	arr[2*i+1]=Out[0].c[2*m+1].col[i];
      }
      float t0=tmax/(npoc-1);
      float norm=tmax/(2*M_PI);
      tmp=four_direct(arr,npoc,1,t0*f_ws,t0*f_we,t0*f_dw);

      for(i=0;i<npof;i++){
	arr[i]=     tmp[2*i  ]*norm;
	arr[dim+i]=tmp[2*i+1]*norm;
      }
      free(tmp);

      //smoothing if needed
      if(out_w>0){
        eprintf("Smoothing functions...\n");
        Smooth(&Fr,iwin);
        Smooth(&Fi,iwin);
      }

      if(f_out[m]&FOUT_R){ // write real part
	sprintf(str2," %d-Re(%s)",k+2,AnDsymb[m]);
	strcat(str,str2);
	Out[1].UpdateColumn(k,  &Fr);
	k++;
      }
      if(f_out[m]&FOUT_I){ // write imaginary part
	sprintf(str2," %d-Im(%s)",k+2,AnDsymb[m]);
	strcat(str,str2);
	Out[1].UpdateColumn(k,  &Fi);
	k++;
      }
      if(f_out[m]&FOUT_M){ // write modulus
	sprintf(str2," %d-Mod(%s)",k+2,AnDsymb[m]);
	strcat(str,str2);
	Fr.operation(Fi,make_mod);
	Out[1].UpdateColumn(k,  &Fr);
	k++;
      }

    }
    delete [] arr;


    UpdateFile(1,str);
    Out[1].dealloc();
  }

  Out[0].dealloc();

  //for(i=0;i<NCORR;i++){
  //  delete [] nd_r[i];
  //}
  return 1;
}





int  AnDstruc::Process(char *name){
  int res=1;
  if(!(psh->sp_type&PROF_SPLIT)){

    char str[256]="dstruc.dat";

    sprintf(str,"%s%s-dy%%s.dat",psh->out_dir,name);
    eprintf("Writing dynamical structure factor...\n");
    res&=fileout(str);

    eprintf("Done.\n");
  }

  return res;
}




int  AnDstruc::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    int i;


    char str[256]="dstruc.dat";

    sprintf(str,"%s%s-dy%%s.%03d",psh->out_dir,psh->dname,nsp);
    eprintf("Writing dynamical structure factor for split #%d...\n",nsp);
    res&=fileout(str);
    eprintf("Done.\n");

    if(!psh->lastsplit){
      eprintf("ReInit: Dynamical structure factor analyser.\n");

      nisteps=psh->regsteps_next(PR_IONCOORDS);
      nesteps=psh->regsteps_next(PR_ELCCOORDS);
      nsteps=nesteps;
      //nisteps=nesteps;

      if(nd_type==VC_FFT){
	eprintf("DSF: Electron Fourier degree: %d",(int)log2(nsteps)+1);
	nd_n=1<<((int)log2(nsteps)+1);
	eprintf(", %d points\n",nd_n);

	eprintf("DSF: Ion Fourier degree: %d",(int)log2(nisteps)+1);
	nd_ni=1<<((int)log2(nisteps)+1);
	eprintf(", %d points\n",nd_ni);
      }
      else{
	nd_n=nsteps;
	nd_ni=nisteps;
      }
      nd_ni=nd_n; // assume for simplicity, linear interpolation

      // distribution renormalization factors (must be equal!!!)
      //factor_i=sqrt((float)nd_ni/nisteps);
      factor_e=1;  //1/psh->split_length(); //   sqrt((float)nd_n/nesteps);
      factor_i=factor_e;

      int m;
      for(m=0;m<NCORR;m++){

	for(i=0;i<nd_nvect;i++){
	  if(m==0){
	    nd_r[m][i].init(nd_ni,psh->split_length(),CORR_COMPLEX);
	    nd_r[m][i].insert_begin(nisteps,factor_i);
	  }
	  else{
	    nd_r[m][i].init(nd_n,psh->split_length(),CORR_COMPLEX);
	    // only for ee and ii
	    if(m<2 || m>3)nd_r[m][i].insert_begin(nsteps,factor_e);
	  }
	}
      }
      cur_isr=cur_isi=0.;
    }
  }
  return res;
}













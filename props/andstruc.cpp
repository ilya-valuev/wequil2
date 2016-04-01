# include <stdio.h>
# include <math.h>

# include "common.h" 
# include "rparam.h"

# ifdef UNIX
# include "../plasma.h"
# include "andstruc.h"
# include "../interact.h"
# include "../plasmar.h"
# include "statist.h"
# else
# include "../plasma.h"
# include "andstruc.h"
# include "../interact.h"
# include "../plasmar.h"
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
  double lsq=n;
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
	f_out[i]|=FOUT_R;
	f_outn++;
      }
      sprintf(str,"%s%s",AnDsymb[i],"i");
      if(strstr(tmpstr,str)){
	f_out[i]|=FOUT_I;
	f_outn++;
      }
      sprintf(str,"%s%s",AnDsymb[i],"m");
      if(strstr(tmpstr,str)){
	f_out[i]|=FOUT_M;
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
    sub_mean=0;
    
    wr_wc=0;
    if(Read_param("DSF write raw square: %s",tmpstr)){
      if(strstr(tmpstr,"y"))
        wr_wc=1;
    }

    Set_stop(1);
  }
  Set_stop(0);
  if(Read_param("DSF subtract mean: %s",tmpstr)){
    if(strstr(tmpstr,"y"))
      sub_mean=1;
  }
  Set_stop(1);

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
  //factor_i=sqrt((double)nd_ni/nisteps);
  factor_e=1; //1./psh->split_length();  //1; //sqrt((double)nd_n/nesteps);
  factor_i=factor_e;

  Set_stop(0);

  int i,l;
  double kval;
  if(Read_param("DSF Free K: %lf %d", &kval, &nd_nvect)==2){
    eprintf("DSF: Warning: using not allowed k-directions! \n");
    nd_nksq=1;
    kind= new int[nd_nksq];
    if(!kind)serror("DSF: MAE\n");

    Vector_3 k;
    for(i=0;i<nd_nvect;i++){
      k=Vector_3(random2(),random2(),random2());
      k.normalize();
      k*=kval*2*M_PI/psh->gasp->L;
      lkvect.Append(k);
    }
    kind[0]=0;
  }
  else{
    Set_stop(1);
    int ksq;

    Read_param("DSF K square:>",tmpstr);
    klist=new cList(tmpstr);
    if(!klist || !klist->valid()){
       msg_error("Can't read valid K vector list (%s)!\n",tmpstr);
       return 0;
    }
    nd_nksq=klist->n;
    kind = new int[nd_nksq];
    if(!kind){
      msg_error("DSF: MAE kind!\n");
      return 0;
    }

    klist->rewind();

    for(i=0;i<nd_nksq;i++){
      kind[i]=nd_nvect;

      ksq=klist->step();
      int nk=splitsq(ksq,3,arr);
      int loc_nv=0;
      Vector_3 k;

      for(l=0;l<nk;l++){ //all possible vectors
        int i;
        for(i=0;i<3;i++){
          k[i]=(double)arr[MAXLEV*l+i];
        }
        k*=2*M_PI/psh->gasp->L;

        // determination of nvect
        int e;
        int m=0;
        while(k[m]==0. && m<3)m++;

        for(e=0;e<4;e++){ // averaging through 8 possible directions

         //nd_kvect[nd_nvect++]=k;
         lkvect.Append(k);
         loc_nv++;
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
      eprintf("K square %d: %d components, %d vectors\n",ksq,nk,loc_nv);
      nd_nvect=lkvect.Count();
    }
  }

  eprintf("Total: %d vectors\n",nd_nvect);
  if(nd_nvect==0)serror("!!!!!!???????\n");

  Set_stop(1);
  single_file=0;
  Read_param("DSF K output:>",tmpstr);
  if(strstr(tmpstr,"one"))single_file=1;

  if(single_file && wr_tc){
    eprintf("DSF: Cannot create single file with time correlations, using multiple files!\n");
    single_file=0;
  }
  if(single_file && f_outn>1){
    eprintf("DSF: Cannot create single file with number of output columns per K\n"
            "     greater than 1, using multiple files!\n");
    single_file=0;
  }
  if(single_file && wr_wc){
    eprintf("DSF: Cannot create single file with raw square, using multiple files!\n");
    single_file=0;
  }

  int m;
  for(m=0;m<NCORR;m++){
   nd_r[m] = new Correlation[nd_nvect];
   if(!nd_r[m])serror("DSF: MAE\n");


   for(i=0;i<nd_nvect;i++){
    if(m==0){
      nd_r[m][i].init(nd_ni,psh->split_length(),CORR_COMPLEX);
      nd_r[m][i].sub_mean=sub_mean;
      nd_r[m][i].insert_begin(nisteps,factor_i);
    }
    else{
      nd_r[m][i].init(nd_n,psh->split_length(),CORR_COMPLEX);
      nd_r[m][i].sub_mean=sub_mean;
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
  if(klist)delete klist;
  if(kind)delete kind;
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
    //psh->anx=psh->gasp->xx;
    psh->xx_used=1;
  }


  int j,l;
  Vector_3 r,k;
  double sr=0,si=0;

  //for(l=0;l<nd_nvect;l++){ // assembling data for  all possible vectors
    //k=nd_kvect[l];
  l=0;
  lkvect.Rewind();
  while(!lkvect.Done()){
    k=lkvect.GetCur();
    lkvect.Next();

    if(psh->status&PR_IONCOORDS){
      sr=si=0.;
      for(j=0;j<nion;j++){ // averaging through nn2 particles
	r=psh->anx[j];

	sr+=cos(r*k);
	si+=sin(r*k);
      }
      sr/=nion;
      si/=nion;
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
    sr/=nelec;
    si/=nelec;
    nd_r[1][l].insert_next(sr,si);
    nd_r[4][l].insert_next(cur_isr-sr,cur_isi-si);
    l++;
  }

  return 1;
}



double make_mod(double x,double y){
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

  // adding outputs in advance
  AddOutput(0.,1.,10,2);
  AddOutput(0.,1.,10,2);
  AddOutput(0.,1.,10,2);

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


  int j,m,k;

  for(k=0;k<nd_nksq;k++){ // loop over k values
    int i_start=kind[k],i_end;
    if(k<nd_nksq-1)i_end=kind[k+1]-1;
    else i_end=nd_nvect-1;
    int n_av=i_end-i_start+1;
    if(n_av==0)continue;

    for(m=0;m<NCORR;m++){ // averaging over vectors
     int n=nd_r[m][i_start].n;

     for(j=0;j<n;j++){
      for(i=i_start+1;i<=i_end;i++){
        nd_r[m][i_start].arr[j]+=nd_r[m][i].arr[j]; // real part
        nd_r[m][i_start].arri[j]+=nd_r[m][i].arri[j]; // imaginary part
        if(nd_type==VC_FFT){
          nd_r[m][i_start].arr[n+j]+=nd_r[m][i].arr[n+j]; // fft^2
          // only for cross-correlations
          if(m>=2)nd_r[m][i_start].arri[n+j]+=nd_r[m][i].arri[n+j]; // imaginary of fft1 mal fft2*
        }
      }
      nd_r[m][i_start].arr[j]/=n_av;  //*n;
      nd_r[m][i_start].arri[j]/=n_av;

      if(nd_type==VC_FFT){
        nd_r[m][i_start].arr[n+j]/=n_av;
        if(m>=2)nd_r[m][i_start].arri[n+j]/=n_av;
      }
     }
    }
  }

  double tmax=nd_r[1][0].maxtime; // tmax and n from <ee>
  double dt=fmax(0.02,tmax/nd_r[1][0].n);

  int npo=(int)(tmax/dt)+1;
  if(npo<2)npo=2;

  TableFunction Sr, Si; // full correlations S=<ii>+<ee>-<ei>-<ie>
  TableFunction Wr, Wi; // full correlations (direct FFT result) S=<ii>+<ee>-<ei>-<ie>



  char str[500]="1-w", str2[250];
  char str_low[500]="1-w";
  k=0; // column number

  double mul=2*M_PI/psh->gasp->L;
  mul*=mul;

  for(i=0;i<nd_nksq;i++){
    if(!single_file){
      k=0;
      strcpy(str,"1-w");
      strcpy(str_low,"1-w");
    }
    
    double k2=(lkvect.Get(kind[i])).norm2()/(mul);  // k square
    if(!ExistsOutput(0))
      serror("Output is not registered!\n");  //AddOutput(0,tmax,npo,12);
    else Out[0].init(0,tmax,npo,12);

    char tst[50]="dt";
    if(wr_tc){
      if(single_file){
        if(i==0){
          tst[0]=tch;
          sprintf(filename,fileform,tst);
          eprintf("Updating time correlations file...\n");
          AssociateOutput(filename,0,psh->update);
        }
        else AssociateOutput(filename,0,0); // to prevent loading
      }
      else{  // multiple files
        sprintf(tst,"%ct_k%.0f",tch,k2);
        sprintf(filename,fileform,tst);
        eprintf("Updating time correlations file for vector %d...\n",i+1);
        AssociateOutput(filename,0,psh->update);
      }
    }
   
    Sr<< nd_r[0][kind[i]].func_r();
    Sr+= nd_r[1][kind[i]].func_r();
    Sr-= nd_r[2][kind[i]].func_r();
    Sr-= nd_r[3][kind[i]].func_r();

    Si<< nd_r[0][kind[i]].func_i();
    Si+= nd_r[1][kind[i]].func_i();
    Si-= nd_r[2][kind[i]].func_i();
    Si-= nd_r[3][kind[i]].func_i();


    Out[0].UpdateColumn(0,&nd_r[0][kind[i]].func_r());  // Re(Sii)
    Out[0].UpdateColumn(1,&nd_r[0][kind[i]].func_i());  // Im(Sii)
    Out[0].UpdateColumn(2,&nd_r[1][kind[i]].func_r());  // Re(See)
    Out[0].UpdateColumn(3,&nd_r[1][kind[i]].func_i());  // Im(See)
    Out[0].UpdateColumn(4,&nd_r[2][kind[i]].func_r());  // Re(Sei)
    Out[0].UpdateColumn(5,&nd_r[2][kind[i]].func_i());  // Im(Sei)
    Out[0].UpdateColumn(6,&nd_r[3][kind[i]].func_r());  // Re(Sie)
    Out[0].UpdateColumn(7,&nd_r[3][kind[i]].func_i());  // Im(Sie)

    Out[0].UpdateColumn(8,&Sr);  // Re(S)
    Out[0].UpdateColumn(9,&Si);  // Im(S)

    Out[0].UpdateColumn(10,&nd_r[4][kind[i]].func_r());  // Re(S), direct sum
    Out[0].UpdateColumn(11,&nd_r[4][kind[i]].func_i());  // Im(S)


    if(wr_tc){
      if(!single_file || i==0){
        UpdateFile(0,"1-t 2-Re(Sii) 3-Im(Sii) 4-Re(See) 5-Im(See) 6-Re(Sie) "
                     "7-Im(Sie) 8-Re(Sei) 9-Im(Sei) 10-Re(S) 11-Im(S) 12-Re(Sdir) 13-Im(Sdir)");
      }
    }
    double tmax;
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
    if(wr_wc){
      //double fw_dw=2*M_PI/nd_r[1][0].maxtime;
      //double fw_we=fw_dw*nd_r[1][0].n;

      //Wr.xscale(0,fw_we);
      //Wi.xscale(0,fw_we);
      // determining smoothing window
      int iwin=1;
      if(out_w>0){
        iwin=(int)(out_w/f_dw+1);
      }

      for(int j=0;j<5;j++){
        //nd_r[j][kind[i]].func_fsqr().xscale(0,fw_we);
        //nd_r[j][kind[i]].func_fsqi().xscale(0,fw_we);
      
        if(iwin>0){ //smoothing if needed
          Smooth(&nd_r[j][kind[i]].func_fsqr(),iwin);
          Smooth(&nd_r[j][kind[i]].func_fsqi(),iwin);
        }
      }
      Wr<< nd_r[0][kind[i]].func_fsqr();
      Wr+= nd_r[1][kind[i]].func_fsqr();
      Wr-= nd_r[2][kind[i]].func_fsqr();
      Wr-= nd_r[3][kind[i]].func_fsqr();

      Wi<< nd_r[0][kind[i]].func_fsqi();
      Wi+= nd_r[1][kind[i]].func_fsqi();
      Wi-= nd_r[2][kind[i]].func_fsqi();
      Wi-= nd_r[3][kind[i]].func_fsqi();
      if(iwin>0){
        Smooth(&Wr,iwin);
        Smooth(&Wi,iwin);
      }

      if(single_file){
        if(i==0){
          if(!ExistsOutput(2))
            serror("Output is not registered!\n"); //AddOutput(0,fw_we,nd_r[1][0].n,12*nd_nksq);
          else Out[2].init(f_ws,f_we,npof,12*nd_nksq);

          tst[0]=tch;
          tst[1]='w';
          sprintf(filename,fileform,tst);
          eprintf("Updating DSF w file...\n");
          AssociateOutput(filename,2,psh->update);
        }
      }
      else{  // multiple files
        if(!ExistsOutput(2))
          serror("Output is not registered!\n"); //AddOutput(0,fw_we,nd_r[1][0].n,12);
        else Out[2].init(f_ws,f_we,npof,12);
        sprintf(tst,"%cw_k%.0f",tch,k2);
        sprintf(filename,fileform,tst);
        eprintf("Updating DSF w file for vector %d...\n",i+1);
        AssociateOutput(filename,2,psh->update);
      }
      

      

      Out[2].UpdateColumn(0,&nd_r[0][kind[i]].func_fsqr());  // Re(Sii)
      Out[2].UpdateColumn(1,&nd_r[0][kind[i]].func_fsqi());  // Im(Sii)
      Out[2].UpdateColumn(2,&nd_r[1][kind[i]].func_fsqr());  // Re(See)
      Out[2].UpdateColumn(3,&nd_r[1][kind[i]].func_fsqi());  // Im(See)
      Out[2].UpdateColumn(4,&nd_r[2][kind[i]].func_fsqr());  // Re(Sei)
      Out[2].UpdateColumn(5,&nd_r[2][kind[i]].func_fsqi());  // Im(Sei)
      Out[2].UpdateColumn(6,&nd_r[3][kind[i]].func_fsqr());  // Re(Sie)
      Out[2].UpdateColumn(7,&nd_r[3][kind[i]].func_fsqi());  // Im(Sie)

      Out[2].UpdateColumn(8,&Wr);  // Re(S)
      Out[2].UpdateColumn(9,&Wi);  // Im(S)

      Out[2].UpdateColumn(10,&nd_r[4][kind[i]].func_fsqr());  // Re(S), direct sum
      Out[2].UpdateColumn(11,&nd_r[4][kind[i]].func_fsqi());  // Im(S)

      char strw[200]="1-t 2-Re(SWii) 3-Im(SWii) 4-Re(SWee) 5-Im(SWee) 6-Re(SWie) "
                     "7-Im(SWie) 8-Re(SWei) 9-Im(SWei) 10-Re(SW) 11-Im(SW) 12-Re(SWdir) 13-Im(SWdir)";
      if(single_file){
        if(i==nd_nksq-1){
          strcat(strw,"\n\n#");
          UpdateFile(2,strw);
          Out[2].dealloc();
        }
      }
      else{  // multiple files
        UpdateFile(2,strw);
        Out[2].dealloc();
      }
    }
    
    if(wr_fc){
      

      if(single_file){
        if(i==0){
          if(!ExistsOutput(1))
            serror("Output is not registered!\n"); //AddOutput(f_ws,f_we,npof,f_outn*nd_nksq);
          else Out[1].init(f_ws,f_we,npof,f_outn*nd_nksq);

          tst[0]=tch;
          tst[1]='f';
          sprintf(filename,fileform,tst);
          eprintf("Updating DSF file...\n");
          if(!wr_tc)AssociateOutput(filename,1,psh->update);
          else AssociateOutput(filename,1,0);
        }
      }
      else{  // multiple files
        if(!ExistsOutput(1))
          serror("Output is not registered!\n"); //AddOutput(f_ws,f_we,npof,f_outn);
        else Out[1].init(f_ws,f_we,npof,f_outn);
        sprintf(tst,"%cf_k%.0f",tch,k2);
        sprintf(filename,fileform,tst);
        eprintf("Updating DSF file for vector %d...\n",i+1);
        if(!wr_tc)AssociateOutput(filename,1,psh->update);
        else AssociateOutput(filename,1,0);
      }

      int dim=fmax(npoc,npof);
      double *arr= new double[2*dim], *tmp;
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
        iwin=(int)(out_w/f_dw+1);
      }

      eprintf("Frequency transform...\n");

      for(m=0;m<6;m++){
        int i;
        if(f_out[m]==FOUT_NONE)continue;

        eprintf("Columns %d and %d...\n",2*m+2, 2*m+3);
        double sumr=0, sumi=0;
        for(i=0;i<npoc;i++){
          arr[2*i  ]=Out[0].c[2*m  ].col[i];
          sumr+=arr[2*i  ];
          arr[2*i+1]=Out[0].c[2*m+1].col[i];
          sumi+=arr[2*i+1];
        }
        sumr/=npoc;
        sumi/=npoc;
        if(sub_mean){
          for(i=0;i<npoc;i++){
            arr[2*i  ]-=sumr;
            arr[2*i+1]-=sumi;
          }
        }

        double t0=tmax/(npoc-1);
        double norm=tmax/(2*M_PI);
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
          if(single_file){
            sprintf(str2," %d-%.0f(%e 1/RDebye)",k+2,k2,sqrt(k2*mul)*RDebye);
            strcat(str_low,str2);
          }
          Out[1].UpdateColumn(k,  &Fr);
          /*if(single_file)*/
          k++;
        }
        if(f_out[m]&FOUT_I){ // write imaginary part
          sprintf(str2," %d-Im(%s)",k+2,AnDsymb[m]);
          strcat(str,str2);
          if(single_file){
            sprintf(str2," %d-%.0f(%e 1/RDebye)",k+2,k2,sqrt(k2*mul)*RDebye);
            strcat(str_low,str2);
          }
          Out[1].UpdateColumn(k,  &Fi);
          /*if(single_file)*/
          k++;
        }
        if(f_out[m]&FOUT_M){ // write modulus
          sprintf(str2," %d-Mod(%s)",k+2,AnDsymb[m]);
          strcat(str,str2);
          if(single_file){
            sprintf(str2," %d-%.0f(%e 1/RDebye)",k+2,k2,sqrt(k2*mul)*RDebye);
            strcat(str_low,str2);
          }
          Fr.operation(Fi,make_mod);
          Out[1].UpdateColumn(k,  &Fr);
          /*if(single_file)*/
          k++;
        }
      }
      delete [] arr;

      if(single_file){
        if(i==nd_nksq-1){
          strcat(str_low,"\n#");
          strcat(str_low,str);
          UpdateFile(1,str_low);
          Out[1].dealloc();
        }
      }
      else{  // multiple files
        UpdateFile(1,str);
        Out[1].dealloc();
      }
    }
    Out[0].dealloc();
  }
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
      //factor_i=sqrt((double)nd_ni/nisteps);
      factor_e=1;  //1/psh->split_length(); //   sqrt((double)nd_n/nesteps);
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













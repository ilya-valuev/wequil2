# include <stdio.h>
# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "anpair.h"
# include "../interact.h"
# include "../plasmar.h"


int AnPair::Init(AnalyseShell *parent,int initype){
  
  if(!Analysator::Init(parent,initype))return 0;


  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&COORD)){
      msg_error("Can't calculate pair distribution function,\n"
		" trajectory file does not contain Coords data.\n");
      return 0;
    }
  }

  Set_stop(0);
  double rr_r0, rr_r1;
  if(Read_param("r-r range: %lf, %lf",&rr_r0,&rr_r1)!=2){
    rr_r0=0.;
    rr_r1=-1.;
  }
  int npo=400;
  Read_param("r-r intervals: %d",&npo);
  rr_coeff=1.;
  Read_param("r-r distance factor: %lf",&rr_coeff);

  char tmpstr[256];
  spin=0;
  if(Read_param("r-r spin: %s",tmpstr)){
    if(!strcmp(tmpstr,"yes"))
     spin=1;
  }


  DRRee.init(rr_r0,(rr_r1>0 ? rr_r1 : psh->gasp->L/2),npo);
  DRRep.init(rr_r0,(rr_r1>0 ? rr_r1 : psh->gasp->L/2),npo); 
  
  if(spin){
    DRReep.init(rr_r0,(rr_r1>0 ? rr_r1 : psh->gasp->L/2),npo);
    DRReeu.init(rr_r0,(rr_r1>0 ? rr_r1 : psh->gasp->L/2),npo); 
  }
  
  if(psh->gasp->non_symm){
    DRRpp.init(0,(rr_r1>0 ? rr_r1 : psh->gasp->L/2),npo);
    eprintf("Init: Non-symmetrical pair distribution analyser.\n");
  }
  else eprintf("Init: Symmetrical pair distribution analyser.\n");

  if(psh->sp_type&AV_SPLIT){
    eprintf("AnPair: warning: split averaging is not implemented!\n");
    psh->sp_type&=~(AV_SPLIT);
  }
  Set_stop(1);
  return 1;
}

int AnPair::Step(){
  if(!(psh->status&PR_IONCOORDS) 
     && !(psh->status&PR_ELCCOORDS))return 1; //nothing to do

  int i,j;
  double R;
  Vector_3 r; 

  int ni=psh->gasp->ni; // only for charge-symmetrical case!
  int npart=psh->gasp->n;
  int nup=psh->gasp->ne/2; 
  

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


  if(psh->status&PR_IONCOORDS){
    for(i=0;i<ni;i++){
      for(j=i+1;j<ni;j++){
	r=psh->anx[j]-psh->anx[i];
	r=psh->gasp->rcell(r);
	R=r.norm();
	if(psh->gasp->non_symm)DRRpp.point(R,1.);
	else DRRee.point(R,1.);
      }
    }
  }	
  
  if(psh->status&PR_ELCCOORDS){
    for(i=ni;i<npart;i++){
      int s1= i-ni<nup ? 1: -1; 
      for(j=i+1;j<npart;j++){
        r=psh->anx[j]-psh->anx[i];
	r=psh->gasp->rcell(r);
	R=r.norm();
	DRRee.point(R,1.);
        if(spin){
          int s2= j-ni<nup ? 1: -1;
          if(s1*s2>0)
            DRReep.point(R,1.);
          else
            DRReeu.point(R,1.);
        }
      }
    }
  }

  if((psh->status&PR_ELCCOORDS) && (psh->status&PR_IONCOORDS)){
    for(i=0;i<ni;i++){
      for(j=ni;j<npart;j++){
	r=psh->anx[j]-psh->anx[i];
	r=psh->gasp->rcell(r);
	R=r.norm();
	DRRep.point(R,1.);
      }
    }
  }

  return 1;
}

double Mfact=1., Mr2=1., Mr1=0., Mr0=0.;

double Mfunc(double x){
  if(fabs(x)<1e-32)return 0.;
  else return Mfact/(Mr2*x*x+Mr1*x+Mr0);
}


int AnPair::WriteDRR(char *file){
  //int i;
  double l=(DRRee.x2-DRRee.x1); 

  int npo=DRRee.n;
  if(npo<2)npo=2;

  /*double na=0.;
  if(DRRee.count)na=l*l*(npo-1)/3/DRRee.count;
  double nb=0.;
  if(DRRep.count)nb=l*l*(npo-1)/3/DRRep.count;
  double nc=0.;
  if(psh->gasp->non_symm && DRRpp.count)nc=l*l*(npo-1)/3/DRRpp.count;

  int nco=4; // number of data columns
  if(psh->gasp->non_symm)nco=6;

  if(!ExistsOutput(0))AddOutput(DRRee.x1,DRRee.x2,npo,nco);
  else Out[0].init(DRRee.x1,DRRee.x2,npo,nco);

  AssociateOutput(file,0,psh->update);

  int ic=0;
  Out[0].UpdateColumn(ic++,&DRRee.distr);
  Mfact=na;
  DRRee.distr*=Mfunc;
  Out[0].UpdateColumn(ic++,&DRRep.distr);

  Mfact=nb;
  DRRep.distr*=Mfunc;
  if(psh->gasp->non_symm){
    Out[0].UpdateColumn(ic++,&DRRpp.distr);
    Mfact=nc;
    DRRpp.distr*=Mfunc;
  }

  Out[0].UpdateColumn(ic++,&DRRee.distr);
  Out[0].UpdateColumn(ic++,&DRRep.distr);
  if(psh->gasp->non_symm){
    Out[0].UpdateColumn(ic++,&DRRpp.distr);
  }*/


  double dl=l/(npo-1);
  double Vtot=psh->gasp->L;
  Vtot=Vtot*Vtot*Vtot;
  Mr2=dl;
  Mr1=0.;
  Mr0=dl*dl*dl/12.;
  double nee=0.;
  if(DRRee.allcount)
    nee=Vtot/(4*M_PI*DRRee.allcount);
  double nii=0.;
  if(psh->gasp->non_symm && DRRpp.allcount)
    nii=Vtot/(4*M_PI*DRRpp.allcount /*psh->sp_nsteps*nii*(nii-1)/2*/);
  double nei=0.;
  if(DRRep.allcount)
    nei=Vtot/(4*M_PI*DRRep.allcount /*psh->sp_nsteps*nei*/);
  double neep=0.;
  if(spin && DRReep.allcount)
    neep=Vtot/(4*M_PI*DRReep.allcount);
  double neeu=0.;
  if(spin && DRReeu.allcount)
    neeu=Vtot/(4*M_PI*DRReeu.allcount);
  

  int nco=4; // number of data columns
  if(psh->gasp->non_symm)nco=6;
  if(spin)
    nco+=2;

  if(!ExistsOutput(0))AddOutput(DRRee.x1,DRRee.x2,npo,nco);
  else Out[0].init(DRRee.x1,DRRee.x2,npo,nco);
  Out[0].SetArgCoeff(rr_coeff);

  AssociateOutput(file,0,psh->update);

  int ic=0;
  Out[0].UpdateColumn(ic++,&DRRee.distr);
  Mfact=nee;
  DRRee.distr*=Mfunc;
  Out[0].UpdateColumn(ic++,&DRRep.distr);

  Mfact=nei;
  DRRep.distr*=Mfunc;
  if(psh->gasp->non_symm){
    Out[0].UpdateColumn(ic++,&DRRpp.distr);
    Mfact=nii;
    DRRpp.distr*=Mfunc;
  }

  Out[0].UpdateColumn(ic++,&DRRee.distr);
  Out[0].UpdateColumn(ic++,&DRRep.distr);
  if(psh->gasp->non_symm){
    Out[0].UpdateColumn(ic++,&DRRpp.distr);
  }

  if(spin){
    Mfact=neep;
    DRReep.distr*=Mfunc;
    Mfact=neeu;
    DRReeu.distr*=Mfunc;
    Out[0].UpdateColumn(ic++,&DRReep.distr);
    Out[0].UpdateColumn(ic++,&DRReeu.distr);
  }

  char header[1000]="", tmp[400];  
  int g=1;
  if(rr_coeff==1.)
    sprintf(tmp,"%d-r/lLandau  %d-F2ee(r) %d-F2ei(r)",g,g+1,g+2);
  else
    sprintf(tmp,"%d-r*d_fact/lLandau  %d-F2ee(r) %d-F2ei(r)",g,g+1,g+2);
  g+=3;
  strcat(header,tmp);
  if(psh->gasp->non_symm){
    sprintf(tmp," %d-F2ii(r) %d-Fee(r)",g,g+1);
    g+=2;
    strcat(header,tmp);
  }
  sprintf(tmp," %d-Fei(r) %d-Fii(r)",g,g+1);
  g+=2;
  strcat(header,tmp);
  
  if(spin){
    sprintf(tmp," %d-Fee_p(r) %d-Fee_u(r)",g,g+1);
    g+=2;
    strcat(header,tmp);
  }
  
  UpdateFile(0,header);
  Out[0].dealloc();


# if 0

 double a,b,c,x,dx,aa,bb,cc;
 dx=(DRRee.x2-DRRee.x1)/399;
 FILE *f1=Err_fopen(file,"wt");
 

 for(x=DRRee.x1,i=0;i<400;x+=dx,i++){
   a=DRRee(x);
   b=DRRep(x);
   if(fabs(x)>1e-10){
    aa=a*na/(x*x);
    bb=b*nb/(x*x);
   }
   else aa=bb=0.;

   if(!psh->gasp->non_symm){
     fprintf(f1,"%f %f %f %f %f \n",x,a*na,b*nb,aa,bb);
   }
   else{
     c=DRRpp(x);
     if(fabs(x)>1e-10)cc=c*nc/(x*x);
     else cc=0.;
     fprintf(f1,"%f %f %f %f %f %f %f \n",x,
     	     a*na,c*nc,b*nb,aa,cc,bb);
   }

 }
 fclose(f1);

# endif

 return 1;
}


int  AnPair::Process(char *name){
  if(!(psh->sp_type&PROF_SPLIT)){
    eprintf("Writing pair distributions...\n");
    char str[256]="statstruc.dat";
    strcpy(str,psh->out_dir);
    strcat(str,name);
    strcat(str,"-rr.dat");
    int res= WriteDRR(str);
    eprintf("Done.\n");
    return res;
  }
  return 1;
}


int  AnPair::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    eprintf("Writing pair distributions for split #%d...\n",nsp);
    char str[256];

    sprintf(str,"%s%s-rr.%03d",psh->out_dir,psh->dname,nsp);
    res=WriteDRR(str);

    if(psh->gasp->non_symm)DRRpp.clear();
    DRRee.clear();
    DRRep.clear();
  }
  return res;
}











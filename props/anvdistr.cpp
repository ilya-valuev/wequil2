# include <stdio.h>
# include "common.h"
# include "rparam.h"
# include "../plasma.h"
# include "anvdistr.h"
# include "../interact.h"
# include "../plasmar.h"


int AnVdistr::Init(AnalyseShell *parent,int initype){

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

  char tmpstr[250];
  double ve0=0.,ve1=2.;
  double vp0=0., vp1=2./sqrt(psh->gasp->mass);
  int gvp=400, gve=400;
  Read_param("Ve distribution range: %lf, %lf",&ve0,&ve1);
  Read_param("Ve distribution grid: %d",&gve);
  Read_param("Vi distribution range: %lf, %lf",&vp0,&vp1);
  Read_param("Vi distribution grid: %d",&gvp);

  av_coord=0;
  if(Read_param("Average coord distributions: %s",tmpstr)){
    if(strstr(tmpstr,"y"))av_coord=1;
  }


  max_comp=0;
  if(Read_param("Compare to Maxwell: %s",tmpstr)){
    if(strstr(tmpstr,"y"))max_comp=1;
  }

  distr_base=BASE_ABS;
  if(Read_param("Distribution base: %s",tmpstr)){
    if(strstr(tmpstr,"abs"))distr_base=BASE_ABS;
    else if(strstr(tmpstr,"common"))distr_base=BASE_COMMON;
    else if(strstr(tmpstr,"separate"))distr_base=BASE_SEPARATE;
  }


  vk=1./(psh->gasp->Wpe*psh->gasp->L);
  Te=Ti=1.;

  int i;
  for(i=0;i<3;i++){
    DVe[i].init(ve0,ve1,gve);
    DVp[i].init(vp0,vp1,gvp);
  }
  DVe[3].init(0.,ve1,gve);
  DVp[3].init(0.,vp1,gvp);

  
  eprintf("Init: velocity distribution analyser.\n");

  if(psh->sp_type&AV_SPLIT){
    eprintf("AnVdistr: warning: split averaging is not implemented!\n");
  }
  Set_stop(1);
  return 1;
}

int AnVdistr::Step(){
  if(!(psh->status&PR_IONVEL)
     && !(psh->status&PR_ELCVEL))return 1; //nothing to do


  int i,j;
  int nion=psh->gasp->ni;
  int npart=psh->gasp->n;

  int tmp=psh->gasp->one_center;
  Vector_3 vcm(0,0,0);
  if(distr_base==BASE_COMMON){
    psh->gasp->one_center=1;
    vcm=psh->gasp->getCMVel();
  }
  else if(distr_base==BASE_SEPARATE){
    psh->gasp->one_center=0;
    vcm=psh->gasp->getCMVel(0); // ion center-of-mass velocity
  }

  Vector_3 v;
  for(i=0;i<nion;i++){
    v=(psh->anv[i]-vcm)/vk;
    for(j=0;j<3;j++){
      DVp[j].point(v[j],1.);
    }
    DVp[3].point(v.norm(),1.);
  }

  if(distr_base==BASE_SEPARATE){
    vcm=psh->gasp->getCMVel(1); // electron center-of-mass velocity
  } // otherwise vcm is not changed

  for(i=nion;i<npart;i++){
    v=(psh->anv[i]-vcm)/vk;
    for(j=0;j<3;j++){
      DVe[j].point(v[j],1.);
    }
    if(DVe[3].point(v.norm(),1.)<0){
      printf("jojo:\n");
    }
  }

  if(max_comp){
    psh->gasp->getT();
    Te=psh->gasp->Te;
    Ti=psh->gasp->Ti;
    STe.next();
    STi.next();
  }

  // restoring temperature measurement parameter
  psh->gasp->one_center=tmp;
  return 1;
}



double maxwell_vx(double m, double kT,double v){
  return sqrt(m/(2*M_PI*kT))*exp(-m*v*v/(2*kT));
}

double maxwell_v2(double m, double kT,double v){
  return 2*m*v*v*maxwell_vx(m,kT,v)/kT;
}


double M_mass=1.;
double M_T=1.;
double M_vk=1;

double Mv(double x){
  return maxwell_vx(M_mass,M_T,x*M_vk)*M_vk;
}

double Mv2(double x){
  return maxwell_v2(M_mass,M_T,x*M_vk)*M_vk;
}



int AnVdistr::write_file(char *fname, Distribution *DV,double mass, double T){

  int npo=DV[0].n;
  int nco=2;
  if(!av_coord)nco=4;
  if(max_comp)nco+=2;
  nco+=2;

  if(!ExistsOutput(0))AddOutput(DV[0].x1,DV[0].x2,npo,nco);
  else Out[0].init(DV[0].x1,DV[0].x2,npo,nco);

  AssociateOutput(fname,0,psh->update);

  int i;
  for(i=0;i<4;i++){
    double k=(double)DV[i].count/DV[i].allcount;
    //printf("i= %d, knorm=%f\n",i,k);
    DV[i].normalize(k);
  }


  int ind=0;
  char tmp[200],header[200]="";
  sprintf(tmp,"#1-v(L*Wpe) ");
  strcat(header,tmp);
  if(!av_coord){
    sprintf(tmp,"%d-vx %d-vy %d-vz ",ind+2,ind+3,ind+4);
    Out[0].UpdateColumn(ind++,&DV[0].distr);
    Out[0].UpdateColumn(ind++,&DV[1].distr);
    Out[0].UpdateColumn(ind++,&DV[2].distr);
  }
  else{
    sprintf(tmp,"%d-<vi> ",ind+2);
    DV[0].distr+=DV[1].distr;
    DV[0].distr+=DV[2].distr;
    DV[0].distr*=1./3;
    Out[0].UpdateColumn(ind++,&DV[0].distr);
  }

  strcat(header,tmp);

  sprintf(tmp,"%d-v^2",ind+2);
  strcat(header,tmp);
  Out[0].UpdateColumn(ind++,&DV[3].distr);

  if(max_comp){
    sprintf(tmp," %d-vx(max) %d-v^2(max)",ind+2,ind+3);
    strcat(header,tmp);

    M_mass=mass;
    M_T=T;
    M_vk=vk;

    Out[0].UpdateColumn(ind++,Mv);
    Out[0].UpdateColumn(ind++,Mv2);
  }

  if(1){
    sprintf(tmp," %d-vx(max1) %d-v^2(max1)",ind+2,ind+3);
    strcat(header,tmp);

    M_mass=mass;
    M_T=DV[3].av2()*vk*vk*mass/3;
    M_vk=vk;

    Out[0].UpdateColumn(ind++,Mv);
    Out[0].UpdateColumn(ind++,Mv2);
  }


  UpdateFile(0,header);
  Out[0].dealloc();

# if 0


  for(x=DV[0].x1,i=0;i<DV[0].n;i++,x+=dx){
    fprintf(f1,"%e ",x);

    if(av_coord){
      double y=0.;
      for(j=0;j<3;j++)y+=DV[j](x);
      y/=3;
      fprintf(f1,"%e ",y);
    }
    else{
      for(j=0;j<3;j++)fprintf(f1,"%e ",DV[j](x));
    }

    if(x<0)fprintf(f1,"0.");
    else fprintf(f1,"%e",DV[3](x));

    if(max_comp){
      double mvx=maxwell_vx(mass,T,x*vk)*vk;
      double mv2=maxwell_v2(mass,T,x*vk)*vk;
      fprintf(f1," %e %e",mvx,mv2);
    }
    fprintf(f1,"\n");
  }
  fclose(f1);

# endif

  return 1;
}


  
 
  

int  AnVdistr::Process(char *name){
  if(!(psh->sp_type&PROF_SPLIT)){
    char str[256]="veldistr.dat";
    eprintf("Writing velocity distributions.\n");
       
    if(max_comp){
      Te=STe.av();
      Ti=STi.av();
    }

    int res=1;
    sprintf(str,"%s%s-ve.dat",psh->out_dir,name);
    res&=write_file(str,DVe,1.,Te);

    sprintf(str,"%s%s-vi.dat",psh->out_dir,name);
    res&=write_file(str,DVp,psh->gasp->mass,Ti);
    eprintf("Done.\n");
    return res;
  }
  return 1;
}


int  AnVdistr::ProcessSplit(){
  int res=1;
  if(psh->sp_type&PROF_SPLIT){
    int nsp=psh->sp_num;
    eprintf("Writing velocity distributions for split #%d...\n",nsp);

    char str[256]="veldistr.dat";
       
    if(max_comp){
      Te=STe.av();
      Ti=STi.av();
    }

 
    sprintf(str,"%s%s-ve.%03d",psh->out_dir,psh->dname,nsp);
    res&=write_file(str,DVe,1.,Te);

    sprintf(str,"%s%s-vi.%03d",psh->out_dir,psh->dname,nsp);
    res&=write_file(str,DVp,psh->gasp->mass,Ti);

    int i;
    for(i=0;i<4;i++){
      DVe[i].clear();
      DVp[i].clear();
    }
    if(max_comp){
      STe.clear();
      STi.clear();
    }
  
  }
  return res;
}






















































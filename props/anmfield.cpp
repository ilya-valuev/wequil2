# include <stdio.h>
# include "anmfield.h"
# include "rparam.h"



int AnMField::Init(AnalyseShell *parent,int initype){

  if(!Analysator::Init(parent,initype))return 0;

  if(!initype){
    if(!psh->rec)return 0;
    if(!(psh->rec->wtype&COORD)){
      msg_error("Can't calculate microfields,\n"
			 "trajectory file does not contain Coords data.\n");
      return 0;
    }
  }

  eprintf("Init: Microfield analyser.\n");
  Set_stop(1);
  int i,j;
  char tmpstr[2000], tmpstr1[2000];
  int ac_n=1;
  int ng[3];

  Read_param("Microfield at: %s",tmpstr);
  if(strstr(tmpstr,"space")){
    //micro_field=2;
    int &np=mc_np;
    //Read_param("Average points: %d %s", &np, tmpstr);
    Read_param("Average points:>",tmpstr1);
    char *grd=strstr(tmpstr1,"grid");
    if(grd){
      if(sscanf(grd+4,"%d %d %d",&ng[0],&ng[1],&ng[2])!=3)
        serror("Invalid grid specification %s\n",tmpstr1);
      np=ng[0]*ng[1]*ng[2];
      strcpy(tmpstr,"grid");
      for(int i=2;i>=0;i--)
        if(ng[i]>1){
          mc_grid=ng[i]; // separator spacing
          break;
        }
    }
    else if(sscanf(tmpstr1,"%d %s", &np, tmpstr)!=2)
      serror("Invalid point specification %s\n",tmpstr1);
    mc_type=ON_SPACE;
    mc_tpx= new Vector_3[np];
  }
  else if(strstr(tmpstr,"ion")){
    mc_np=psh->gasp->ni;
    mc_type=ON_ION;
    mc_tpx=psh->anx+psh->gasp->ne;
    
  }
  else if(strstr(tmpstr,"electron")){
    mc_np=psh->gasp->ne;
    mc_type=ON_ELC;
    mc_tpx=psh->anx;
    
  }
  else serror("Unrecognized microfield calculation type: %s\n",tmpstr);
 
  typedef char *charP;
  mc_tpname= new charP[mc_np];

  if(!mc_tpx || !mc_tpname)serror("Error allocating memory\n");
  for(i=0;i<mc_np;i++){
    mc_tpname[i]=new char[20];
    if(!mc_tpname[i])serror("Error allocating memory\n");
  }
  
  if(mc_type==ON_SPACE){
    if(strstr(tmpstr,"random")){
      mc_rnd=1;
      for(i=0;i<mc_np;i++){
        sprintf(mc_tpname[i],"point%d",i);
        mc_tpx[i]=Vector_3(psh->gasp->L*random2()/2,psh->gasp->L*random2()/2,psh->gasp->L*random2()/2);
      }
    }
    else if(strstr(tmpstr,"coords")){
      Set_position("Test coords:");
      for(i=0;i<mc_np;i++){
        Read_param("$: %lf %lf %lf",&mc_tpx[i][0],&mc_tpx[i][1],&mc_tpx[i][2]);
        strncpy(mc_tpname[i],AcFormat,20);
      }
    }
    else if(strstr(tmpstr,"grid")){
      for(int i=0;i<ng[0];i++)
        for(int j=0;j<ng[1];j++)
          for(int k=0;k<ng[2];k++){
            int ind=i*ng[1]*ng[2]+j*ng[2]+k;
            sprintf(mc_tpname[ind],"grd(%d,%d,%d)",i,j,k);
            mc_tpx[ind]=Vector_3( ng[0]> 1 ? psh->gasp->L*(((double)i)/(ng[0]-1)-0.5) : 0. ,
                                                       ng[1]> 1 ? psh->gasp->L*(((double)j)/(ng[1]-1)-0.5) : 0. ,
                                                        ng[2]> 1 ? psh->gasp->L*(((double)k)/(ng[2]-1)-0.5) : 0. );
          }
          
          
    }
    else serror("Unrecognized specification %s\n",tmpstr);
    
    Set_stop(0);
    mc_ntype=0;
    if(Read_param("Space_particle: %s"),tmpstr){
      if(strstr(tmpstr,"electron"))mc_ntype=0x1;
    }
    Set_stop(1);
      
	    
  }
  else if(mc_type==ON_ION){
    for(i=0;i<mc_np;i++)sprintf(mc_tpname[i],"ion%d",i);
  }
  else if(mc_type==ON_ELC){
    for(i=0;i<mc_np;i++)sprintf(mc_tpname[i],"electron%d",i);
  }

  Read_param("Distribution range: %lf, %lf, %d",&mc_r0, &mc_r1, &mc_nr);
  mc_F0=2.603*pow(psh->gasp->n/(psh->gasp->L*psh->gasp->L*psh->gasp->L),2./3.); // base microfield
  
  Read_param("Print results for points: %s", tmpstr);
  if(mc_plist)delete mc_plist;
  if(strstr(tmpstr,"no"))mc_plist=new /*(cList &)*/cList(0,0);
  else if(strstr(tmpstr,"all"))mc_plist=new /*(cList &)*/cList(0,mc_np);
  else mc_plist=new /*(cList &)*/cList(tmpstr);

  //Read_param("Print average for all points: %s", tmpstr);
  //if(strstr(tmpstr,"y"))mc_prav=1;

  Read_param("Print coord distributions: %s",tmpstr);
  if(strstr(tmpstr,"y"))mc_prcoord=1;

  Read_param("Create files for gnuplot: %s",tmpstr);
  if(strstr(tmpstr,"y"))mc_xmgr=1;

  Set_stop(0);
  if(Read_param("Microfield autocorrelation: %s",tmpstr)){
    if(strstr(tmpstr,"y")){
      printf("Microfield: fourier degree: %d",(int)log2(psh->nsteps)+1);
      ac_n=1<<((int)log2(psh->nsteps)+1);
      printf(", %d points\n",ac_n);
      mc_correl=1;     
    }
  }
  mc_cut=-1;
  Read_param("Microfield cutoff distance: %lf",&mc_cut);

  mc_from=0;
  if(Read_param("Microfield from: %s",tmpstr)){
    if(strstr(tmpstr,"ion")) 
      mc_from=ON_ION;
    else if(strstr(tmpstr,"electron"))
      mc_from=ON_ELC;
    else if(!strstr(tmpstr,"all"))
      serror("Unrecognized setting: Microfield from: %s\n",tmpstr);
  }

  
  mc_force=0x1;
  if(Read_param("Force distribution:>",tmpstr)){
    if(strstr(tmpstr,"coulomb"))mc_force|=COULOMB;
    if(strstr(tmpstr,"effect"))mc_force|=EFFECTIVE;
    if(strstr(tmpstr,"all"))mc_force|=ALL;
  }

  mc_anim=0;
  if(Read_param("Microfield animation:>",tmpstr)){
    if(strstr(tmpstr,"yes")){
      mc_anim=1;
    }
  }

  Set_stop(1);

  if(mc_force&EFFECTIVE){
     potspec_t reader;
     reader.read_spec(Get_cur_filename());
     ion_charge=psh->gasp->q;
     mc_potential=reader.potential;
     strncpy(psh->gasp->charpot,reader.charpot,50);


     /*Read_param("Potential: %s",tmpstr);
     if(strstr(tmpstr,"Kelbg"))mc_potential=PotentialKELBG;
     else if(strstr(tmpstr,"Lennard-Johhnes"))mc_potential=PotentialJONES;
     else if(strstr(tmpstr,"Erf"))mc_potential=PotentialERF;
     else if(strstr(tmpstr,"Cutoff")){
       mc_potential=PotentialCUT;
       Read_param("Cutoff value*: %lf",&E_cut);
     }
     else if(strstr(tmpstr,"ln")){
       mc_potential=PotentialLN;
       Read_param("Cutoff value*: %lf",&E_cut);
     }
     else if(strstr(tmpstr,"table")){
       mc_potential=PotentialTAB;
       Read_param("Potential table file: %s",tmpstr);
       char cfgfile[250];
       strcpy(cfgfile,Get_cur_filename());
       Close_param_file();
       ReadPotential(tmpstr);
       Open_param_file(cfgfile);
     }
     else serror("Unknown potential type specified!\n");*/

     reader.calc_lambda(psh->gasp->par_T,1,1,psh->gasp->mass);
     reader.calc_correction(psh->gasp->par_T);
     char str[2000];
     strcpy(str,psh->out_dir);
     //strcat(str,dname);
     strcat(str,"mfeff.pot");
     reader.write_pot(str,psh->gasp->L);
     
     /*int pot_corr=0;
     Read_param("Potential correction: %s",tmpstr);
     if(strstr(tmpstr,"y"))pot_corr=1;
     double Clam_ee=1.,Clam_ep=1.;
     
     double Lambda=0.179*sqrt(psh->gasp->T);
     equal_lambda=1;
  
     Lambda_pauli=0.179*sqrt(psh->gasp->T);
     //non_symm=0;
     Lambda_ee=Lambda_pp=Lambda*Clam_ee;
     if(fabs(Clam_ee-1.)>1e-5)equal_lambda=0;
     Lambda_ep=Lambda*Clam_ep;

     if(pot_corr){
       Correction(IONION,mc_potential,psh->gasp->T);
       Correction(IONELC,mc_potential,psh->gasp->T);
       Correction(ELCELC,mc_potential,psh->gasp->T);
     }*/


  }

  if(mc_force&COULOMB){
    Mf=new Distrvect[mc_np];
    if(!Mf)serror("Microfield_init: memory allocation error.\n");
  }
  if(mc_force&EFFECTIVE){
    Ef=new Distrvect[mc_np];
    //Ef=new Distrvect[mc_np];
    if(!Ef)serror("Microfield_init: memory allocation error.\n");
  }


  for(j=0;j<mc_np;j++){
    for(i=0;i<3;i++){
      if(mc_force&COULOMB)Mf[j][i].init(mc_r0,mc_r1,mc_nr);
      if(mc_force&EFFECTIVE)Ef[j][i].init(mc_r0,mc_r1,mc_nr);
    }
    if(mc_force&COULOMB)Mf[j][3].init(0,mc_r1*sqrt(3.),mc_nr);
    if(mc_force&EFFECTIVE)Ef[j][3].init(0,mc_r1*sqrt(3.),mc_nr);
  }

  if(mc_correl){
    if(mc_force&COULOMB){
      MfCorr=new Correlation[3*mc_np];
      if(!MfCorr)serror("Microfield_init: memory allocation error.\n");
    }
    if(mc_force&EFFECTIVE){
      EfCorr=new Correlation[3*mc_np];
      if(!EfCorr)serror("Microfield_init: memory allocation error.\n");
    }

    for(i=0;i<mc_np;i++){
      for(j=0;j<3;j++){
        if(mc_force&COULOMB){
          MfCorr[3*i+j].init(ac_n,psh->nsteps*psh->dt);
          MfCorr[3*i+j].insert_begin(psh->nsteps);
        }
        if(mc_force&EFFECTIVE){
          EfCorr[3*i+j].init(ac_n,psh->nsteps*psh->dt);
          EfCorr[3*i+j].insert_begin(psh->nsteps);
        }
      }
    }
  }
  return 1;
}



int AnMField::Step(){
  if(!(psh->status&PR_IONCOORDS)
     && !(psh->status&PR_ELCCOORDS))return 1; //nothing to do

  int nion=psh->gasp->ni;
  int nelec=psh->gasp->ne;
  int npart=nion+nelec;

  int i, j;
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

  Vector_3 fe(0,0,0), fc(0,0,0),r;
  double R,f;
  double dEpotent, dQuant;

  
  int type=0;
  if(mc_type==ON_ELC){
    mc_tpx=psh->gasp->xx+nion; //psh->anx+nion;
    type=0x1;
  }
  else if(mc_type==ON_ION){
    mc_tpx=psh->gasp->xx; //psh->anx;
    type=0;
  }
  else // ON_SPACE
    type=mc_ntype;
  
  
  FILE *fa=NULL;
  if(mc_anim){ // animation files
    char str[2000], frm[2000];
    strcpy(str,psh->out_dir);
    //strcat(str,dname);
    strcat(str,"mf%04d.d");
    sprintf(frm,str,mc_anim++);
    fa=fopen(frm,"wt");
    fprintf(fa,"#1-x 2-y 3-z 4-Fcx 5-Fcy 6-Fcz 7-Fex 8-Fey 9-Fez\n");
  }

  for(j=0;j<mc_np;j++){
  
    fe=Vector_3(0,0,0);
    fc=Vector_3(0,0,0);

    int istart=0, iend=npart;
    if(mc_from==ON_ION){
      istart=0;
      iend=nion;
    }
    else if(mc_from==ON_ELC){
      istart=nion;
      iend=npart;
    }
    
    for(i=istart;i<iend;i++){

      if(i>=nion)type|=0x2;// setting ?-electron interaction type 
      
      if(mc_type==ON_ELC && i==j+nion)continue;
      if(mc_type==ON_ION && i==j)continue;
       
      //r=psh->anx[i]-mc_tpx[j];
      r=psh->gasp->xx[i]-mc_tpx[j];
      
      int k;
      for(k=0;k<3;k++){
	      if(r[k]>psh->gasp->L/2)r[k]-=psh->gasp->L;
	      if(r[k]<-psh->gasp->L/2)r[k]+=psh->gasp->L;
      }

    
     
      //r.print();
      R=r.normalize();
      if(mc_cut>0 && R<mc_cut)
        continue;
      if(R>1e-10){
	      if(mc_force&COULOMB){
	        double a=R*R*mc_F0;
	        if(i<nelec)fc+=r/a;
	        else fc-=r/a;
	      }
	      if(mc_force&EFFECTIVE){
	        f=mc_potential(type,R,dEpotent,dQuant);
	        fe+=r*f/mc_F0;
	      }
      }
      else printf("(%d %d)!\n",j,i);
    }

    if(mc_force&COULOMB){
      for(i=0;i<3;i++)Mf[j][i].point(fc[i],1.);
      Mf[j][3].point(fc.norm(),1.);
    }
    if(mc_force&EFFECTIVE){
      for(i=0;i<3;i++)Ef[j][i].point(fe[i],1.);
      Ef[j][3].point(fe.norm(),1.);
    }

    if(mc_correl){
      if(mc_force&COULOMB){
	      for(i=0;i<3;i++)MfCorr[3*j+i].insert_next(fc[i]);
      }
      if(mc_force&EFFECTIVE){
	      for(i=0;i<3;i++)EfCorr[3*j+i].insert_next(fe[i]);
      }
    }

    if(mc_anim){
      fprintf(fa,"%g %g %g %g %g %g %g %g %g\n",mc_tpx[j][0]/psh->gasp->L,mc_tpx[j][1]/psh->gasp->L,mc_tpx[j][2]/psh->gasp->L,fc[0],fc[1],fc[2],fe[0],fe[1],fe[2]);
      if(mc_grid>0 && (j+1)%mc_grid==0)
        fprintf(fa,"\n"); // separator
    }

  }

  if(mc_anim)
    fclose(fa);

  return 1;
}   

# if 0

double Holtsmark(double x){
  if(x<=1e-6)return 0.;

  double dt=0.01*2*M_PI/x;
  double t=dt,sum=0.,ds;
  do{
   ds=exp(-pow(t,3./2.))*sin(x*t);
   sum+=ds;
  }while(fabs(ds)>1e-5);

  return (2/M_PI)*x*sum;
}

# endif

void AnMField::output_mfdistr(char *filename, char *dname, Distrvect *Mf){
  char str[250];
  int i;
  FILE *f1=Err_fopen(filename,"wt"), *f2=NULL;

  if(mc_xmgr){
    strcpy(str,psh->out_dir);
    strcat(str,dname);
    strcat(str,"-mf.gnu");
    f2=Err_fopen(str,"wt");
  }
  
    // creating header
  fprintf(f1,"#// Microfield distribution data\n");
  if(mc_type==ON_SPACE){
    fprintf(f1,"#Test points (coords in cell length L):\n");
    for(i=0;i<mc_np;i++){
      fprintf(f1,"#%s: %f, %f, %f\n",mc_tpname[i],mc_tpx[i][0]/psh->gasp->L,mc_tpx[i][1]/psh->gasp->L,mc_tpx[i][2]/psh->gasp->L);
    }
    fprintf(f1,"\n\n#%13s","   ");
  }
  int j;
  int lfield=13;
  if(mc_prcoord)lfield=13*4;
  char frm[20];
  sprintf(frm,"%%%ds",lfield);
  while((j=mc_plist->step())>=0){
    fprintf(f1,frm,mc_tpname[j]);
  }
  //fprintf(f1,frm,"average");
  //fprintf(f1,frm,"average_c");
  //fprintf(f1,frm,"Holtsmark");
  fprintf(f1,"\n#%s","1-F/F0");

  int ind=2;
  char xyz[]="W(Fx)\0W(Fy)\0W(Fz)\0r";
  for(i=0;i< mc_plist->n + 1;i++){
    char str[50]="";
    if(i!=0)
      sprintf(str,"(point%d)",i);
    if(mc_prcoord)
      for(j=0;j<3;j++)fprintf(f1," %d-%s%s",ind++,xyz+6*j,str);
    //fprintf(f1," %d-%s(p%d)",ind++,"Wr*F^2",i);
  }
  
  
  fprintf(f1," %d-%s",ind++,"4*PI*W(F)*F^2");
  fprintf(f1," %d-%s",ind++,"Wc(F)");
  fprintf(f1," %d-%s",ind++,"4*PI*Nmes(F)*F^2");
  
  //fprintf(f1,"%13s","Wr*F^2");
  fprintf(f1,"\n");
  
  Statistics st[6];
 
 
  double dx, lim=mc_r1*sqrt(3.);
  double x, x0=(mc_prcoord ? mc_r0 : 0);

  /*   No normalization !!!
  for(i=0;i<mc_np;i++){
    for(j=0;j<4;j++)Mf[i][j].normalize(1.);
    //Mf[i][2].normalize_r2(1/(4*M_PI));
  } */
  

  dx=(lim-x0)/mc_nr; // /5;
  
  for(x=x0;x<=lim;x+=dx){

   for(j=0;j<6;j++)st[j].clear();
   for(i=0;i<mc_np;i++){
     for(j=0;j<3;j++){
      st[j].next(Mf[i][j].nrma(x),1);
      st[4].next(Mf[i][j].nrma(x),1);
      //st[4].next(Mf[i][j].nrm_r2(x)+Mf[i][j].nrm_r2(-x),1);
     }
     st[3].next(Mf[i][3].nrma(x),1);
     st[5].next(Mf[i][3](x),1);
   }
  
   fprintf(f1,"%.12g ",x);
   if(mc_xmgr)fprintf(f2,"%.12g ",x);

   mc_plist->rewind();
   while((j=mc_plist->step())>=0){
     if(mc_prcoord){
       for(i=0;i<3;i++){
	 fprintf(f1,"%.12g ",(Mf[j][i])(x));
	 if(mc_xmgr)fprintf(f2,"%.12g ",(Mf[j][i])(x));
       }
       //fprintf(f1,"%12f ",(Mf[2][i])(x)*x*x);
       //if(mc_xmgr)fprintf(f2,"%12f ",(Mf[2][i])(x)*x*x);
     }
     fprintf(f1,"%.12g ",(Mf[j][3])(x));
     if(mc_xmgr)fprintf(f2,"%.12g ",(Mf[j][3])(x));
   }
   if(mc_prcoord){
     for(i=0;i<3;i++){
       fprintf(f1,"%.12g ",st[i].av());
       if(mc_xmgr)fprintf(f2,"%.12g ",(double)st[i].av());
     }
     //fprintf(f1,"%12f ",st[2].av()*x*x);
     //if(mc_xmgr)fprintf(f2,"%12f ",st[2].av()*x*x);
   }
   fprintf(f1,"%.12g ",st[3].av());
   if(mc_xmgr)fprintf(f2,"%.12g ",(double)st[3].av());
   fprintf(f1,"%.12g ",st[4].av());
   if(mc_xmgr)fprintf(f2,"%.12g ",(double)st[4].av());
   fprintf(f1,"%.12g\n",st[5].av());
   if(mc_xmgr)fprintf(f2,"%.12g\n",(double)st[5].av());
   //fprintf(f1,"%12f  \n",Holtsmark(x));
   //if(mc_xmgr)fprintf(f2,"%12f \n",Holtsmark(x));
  }

  if(mc_xmgr)fclose(f2); 
  fclose(f1);
}

void AnMField::output_mfcorr(char *filename,Correlation* MfCorr){
  int i;
  FILE *f1=Err_fopen(filename,"wt");
  printf("Microfield corelations: performing  FFT...\n");
  double *fft_sq=Calculate_forth(MfCorr[0].n,MfCorr,3*mc_np);
  printf("Inverse FFT...\n");
  Calculate_back(fft_sq,MfCorr[0].n);
  for(i=0;i<MfCorr[0].n;i++)fft_sq[i]=fft_sq[2*i];
  TableFunction acf(MfCorr[0].n,fft_sq);
  double tmax=MfCorr[0].maxtime*psh->gasp->Wpe;
  acf.xscale(0.,tmax);
  double dt=0.05,t;
  for(t=0.;t<tmax && t< 150.;t+=dt){
    fprintf(f1,"%f  %f\n",t,acf(t)/acf(0.));
  }
  fclose(f1);
}

int AnMField::Process(char *dname){
  char str[256]="microf.dat";
 
  if(mc_force&COULOMB){
    strcpy(str,psh->out_dir);
    strcat(str,dname);
    if(mc_from==ON_ION)
      strcat(str,"-i");
    else if(mc_from==ON_ELC)
      strcat(str,"-e");
    if     (mc_type==ON_ION)strcat(str,"-mfi.dat");
    else if(mc_type==ON_ELC)strcat(str,"-mfe.dat");
    else strcat(str,"-mfn.dat");
    output_mfdistr(str,dname, Mf);
    
    if(mc_correl){
      strcpy(str,psh->out_dir);
      strcat(str,dname);
      if(mc_from==ON_ION)
        strcat(str,"-i");
      else if(mc_from==ON_ELC)
        strcat(str,"-e");
      if     (mc_type==ON_ION)strcat(str,"-EEi.dat");
      else if(mc_type==ON_ELC)strcat(str,"-EEe.dat");
      else strcat(str,"-EEn.dat");
      output_mfcorr(str,MfCorr);   
    }
  }  
  if(mc_force&EFFECTIVE){
    strcpy(str,psh->out_dir);
    strcat(str,dname);
    if(mc_from==ON_ION)
      strcat(str,"-i");
    else if(mc_from==ON_ELC)
      strcat(str,"-e");
    if     (mc_type==ON_ION)strcat(str,"-efi.dat");
    else if(mc_type==ON_ELC)strcat(str,"-efe.dat");
    else strcat(str,"-efn.dat");
    output_mfdistr(str,dname, Ef);
    
    if(mc_correl){
      strcpy(str,psh->out_dir);
      strcat(str,dname);
      if(mc_from==ON_ION)
        strcat(str,"-i");
      else if(mc_from==ON_ELC)
        strcat(str,"-e");
      if     (mc_type==ON_ION)strcat(str,"-ffi.dat");
      else if(mc_type==ON_ELC)strcat(str,"-ffe.dat");
      else strcat(str,"-ffn.dat");
      output_mfcorr(str,EfCorr);   
    }
  }
  return 1;
}  

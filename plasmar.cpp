# include <stdio.h>
# include <string.h>

# include "common.h"
# include "rparam.h"
# include "plasmar.h"              
# include "plasma.h"



void WriteTypeOut(char *str,int wtype){
  str[0]=0;
  if(wtype&COORD)strcat(str,"coords");
  if(wtype&VEL){
    if(str[0])strcat(str,",");
    strcat(str,"velocities");
  }
  if(wtype&FLOW){
    if(str[0])strcat(str,",");
    strcat(str,"flow");
  }
			
  if(!str[0])strcat(str,"nothing");
}

 
int WriteTypeIn(char *str){
  int type=0;
  if(strstr(str,"coords"))type|=COORD;
  if(strstr(str,"velocities"))type|=VEL;
  if(strstr(str,"flow"))type|=FLOW;
  return type;
}


void PlasmaRecord::write_charhead(){
  Record::write_charhead();
  if(plasmap==NULL)return;
  fprintf(recfp,"\nPlasma parameters:\n");
  plasmap->fileout(recfp);
  fprintf(recfp,"\nRecord parameters:\n");
  fprintf(recfp,"Time step: %f (%f Wpe^-1)\n",plasmap->dt,
	  plasmap->dt*plasmap->Wpe);
  char str[150];
  WriteTypeOut(str,wtype);
  fprintf(recfp,"Write type: %s\n",str);
  fprintf(recfp,"Ion write interval: %ld\n",wions);
  fprintf(recfp,"Electron write interval: %ld\n",welec);
  fprintf(recfp,"N of electron writes per ion write: %ld\n",wenseq);
}





int PlasmaRecord::Init(char *file,int type,Plasma &Gas,long wr_int,
		       long wr_ions,long wr_enseq){
  if(!SetFileName(file))return 0;
  wtype=type;
  
  plasmap=&Gas;

  long wr_idle=wr_ions-wr_enseq;
  if(wr_idle<0){
    eprintf("PlasmaRecord.Init: electron record length (%ld) is"
	    "greater than ions write interval (%ld)!\n"
	    "Setting instanteneous electrons write!\n",wr_enseq,wr_ions);
    wr_idle=0;
    wr_enseq=-1;
  }
  if(wr_int==0)wr_int=-1;
  if(wr_ions==0)wr_ions=-1;



  welec=wr_int;
  wions=wr_ions;
  wenseq=wr_enseq;
  int r;

  if(wtype&COORD){
    if(wions>0){
      r=RegisterFrame("ion_coords",Gas.ni*sizeof(Vector_3),wr_ions);
      if(r<=0)return 0;
      frame[r-1].spec=PR_IONCOORDS;
    }
    if(welec>0){
      r=RegisterFrame("electron_coords",Gas.ne*sizeof(Vector_3),
		      wr_int,0,wr_enseq,wr_idle);
      if(r<=0)return 0;
      frame[r-1].spec=PR_ELCCOORDS;
    }
  }
  
  if(wtype&VEL){
    if(wions>0){
      r=RegisterFrame("ion_velocities",Gas.ni*sizeof(Vector_3),wr_ions);
      if(r<=0)return 0;
      frame[r-1].spec=PR_IONVEL; 
    }
    if(welec>0){
      r=RegisterFrame("electron_velocities",Gas.ne*sizeof(Vector_3),
		      wr_int,0,wr_enseq,wr_idle);
      if(r<=0)return 0;
      frame[r-1].spec=PR_ELCVEL; 
    }
  }

  if(wtype&FLOW){
    r=RegisterFrame("flow",Gas.ne*sizeof(Vector_3),
		      wr_int,0,wr_enseq,wr_idle);
    if(r<=0)return 0;
    frame[r-1].spec=PR_FLOW;
  }

  /*
  // registering reload frames
  if(!wions || !(wtype&COORD)){
     if(RegisterFrame("reload ion_coords",Gas.ni*sizeof(Vector_3),
		     wr_int,0,1,-1)<0)return 0;
  }
  if(!wions || !(wtype&VEL)){
     if(RegisterFrame("reload ion_velocities",Gas.ni*sizeof(Vector_3),
		     wr_int,0,1,-1)<0)return 0;
  }
  if(!welec || !(wtype&VEL)){
     if(RegisterFrame("reload electron_velocities",Gas.ni*sizeof(Vector_3),
		     wr_int,0,1,-1)<0)return 0;
  }
  if(!welec || !(wtype&COORD)){
     if(RegisterFrame("reload electron_velocities",Gas.ni*sizeof(Vector_3),
		     wr_int,0,1,-1)<0)return 0;
  }*/
  
  dstep=1.;
  real_step=(1e-5-1.)*dstep;
  valid=1;
  return 1;
}
  

int PlasmaRecord::Open(char *filename){
  if(Record::Open(filename)<0)return -1;

  Set_scan_end('@');
  Open_param_file(filename);
  Set_stop(1);
  Read_param("Ne: %d",&wne);
  Read_param("Ni: %d",&wni);
  Read_param("Ion mass: %lf\n",&wmass);

  Set_stop(0);
  if(!Read_param("Ion charge: %d\n",&wq)){
    wq=1;
  }
  Set_stop(1);

  Read_param("Potential: %s\n",wpot);
  Read_param("T: %lf",&wT);
  wT*=1e-4;
  Read_param("Gamma: %lf",&wGamma);

  char str[50];
  Set_stop(0);
  wnsymm=1;
  if(Read_param("Symmetric: %s",str)){
    if(strstr(str,"yes"))wnsymm=0;
  }
  wdname[0]=0;
  Read_param("Data name: %s",wdname);
  Set_stop(1);

  Read_param("Write type:>",str);
  wtype=WriteTypeIn(str);


  Read_param("Time step: %lf",&wstp);
  Read_param("Ion write interval: %ld\n",&wions);
  Read_param("Electron write interval: %ld\n",&welec);
  Read_param("N of electron writes per ion write: %ld\n",&wenseq);
  Close_param_file();
  Set_scan_end(0);
  valid=1;

  // setting frame specificators
  int i;
  for(i=0;i<nfr;i++){
    if(strstr(frame[i].name,"ion_velocities"))frame[i].spec=PR_IONVEL;
    else if(strstr(frame[i].name,"electron_velocities"))frame[i].spec=PR_ELCVEL;
    else if(strstr(frame[i].name,"electron_coords"))frame[i].spec=PR_ELCCOORDS;
    else if(strstr(frame[i].name,"ion_coords"))frame[i].spec=PR_IONCOORDS;
    else if(strstr(frame[i].name,"flow"))frame[i].spec=PR_FLOW;
    else frame[i].spec=PR_UNKNOWN;
  }

  dstep=1.;
  real_step=(step+1e-5)/dstep;
  return 0;
}




int PlasmaRecord::Check(char *filename,int type,Plasma &Gas,long wr_int,
		       long wr_ions,long wr_enseq, int num_only){
  if(Open(filename)){
    msg_error("Cannot initialize record: %s\n",filename);
    return 0;
  }

  plasmap=&Gas;

  long wr_idle=wr_ions-wr_enseq;
  if(wr_idle<0){
    eprintf("PlasmaRecord.Check: electron record length (%ld) is"
	    "greater than ions write interval (%ld)!\n"
	    "Setting instanteneous electrons write!\n",wr_enseq,wr_ions);
    //wr_idle=0;
    wr_enseq=-1;
  }
  if(wr_int==0)wr_int=-1;
  if(wr_ions==0)wr_ions=-1;

  // comparing
  eprintf("Comparing:\n"
         "%15s%15s%15s\n","","Init","File");
  int cmp=0;

  eprintf("%15s%15s%15s\n","Potential",Gas.charpot,wpot);

  eprintf("%15s%15d%15d","Ne",Gas.ne,wne);
  if(Gas.ne !=wne){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");

  eprintf("%15s%15d%15d","Ni",Gas.ni,wni);
  if(Gas.ni !=wni){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");

  eprintf("%15s%15.0lf%15d","Qi",Gas.q,wq);
  if(Gas.q !=wq){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");


  eprintf("%15s%15d%15d","Type",type,wtype);
  if(type!=0 && type!=wtype){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");


  eprintf("%15s%15f%15f","Gamma",Gas.par_Gamma,wGamma);
  if(fabs(Gas.par_Gamma-wGamma)>1e-5 && !num_only){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");

  eprintf("%15s%15f%15f","T",Gas.par_T,wT);
  if(fabs(Gas.par_T-wT)>1e-5 && !num_only){
    cmp++;
    eprintf(" !\n");
  }
  else eprintf("\n");



  eprintf("%15s%15f%15f\n","Step",Gas.dt,wstp);

  if(wtype!=0){

    eprintf("%15s%15ld%15ld\n","i-write:",wr_ions,wions);

    eprintf("%15s%15ld%15ld","e-write:",wr_int,welec);
    // cheking proportionality
    double r1=(double)wr_ions/wr_int, r2=(double)wions/welec;
    if(fabs((r1-r2)/r1)>0.01 && !num_only){
      cmp++;
      eprintf(" non-prop\n");
    }
    else eprintf("\n");

    eprintf("%15s%15ld%15ld","e-Nseq:",wr_enseq,wenseq);
    if(wr_enseq!=wenseq && !num_only){
      cmp++;
      eprintf(" !\n");
    }
    else eprintf("\n");
  }


  if(cmp!=0){
    msg_error("Trajectory '%s' is inconsistent"
	      " with Init data.\n",filename);
    return 0;
  }

  valid=1;
  return 1;
}


double PlasmaRecord::AdjustInterval(double dt){

  if(fabs(wstp-dt)>1e-10){
    eprintf("Adjusting calculation step and write interval:\n");
    eprintf("Old step/required step: %f / %f\n",wstp,dt);
  }
  else return wstp;

  int i;
  long *arr= new long[3*nfr];
  if(!arr){
    msg_error("PlasmaRecord::AdjustInterval: MAE!\n");
    return wstp;
  }

  int k=0;
  for(i=0;i<nfr;i++){
    // listing all non-special frames
    if(!strstr(frame[i].name,"reload")){
      // all possible tick intervals
      arr[k++]=frame[i].fq;
      if(frame[i].delay>0)arr[k++]=frame[i].delay;
      if(frame[i].idle>0)arr[k++]=frame[i].idle;
    }
  }
  
  if(k==0){
    delete arr;
    return dt;
  }

  long maxstp=com_divl(k,arr);

  //printf("hehe: %f \n",(wstp*(double)maxstp/dt+0.001));
  long N=(long)(wstp*(double)maxstp/dt+0.1);
  if(N==0){
    eprintf("Step is not adjustable, keeping old.\n");
    delete arr;
    return wstp;
  }
  
  dt=wstp*maxstp/N;
  dstep=(double)((double)maxstp/N+1e-5);
  real_step=step/dstep;

  /*
  for(i=0;i<nfr;i++){
    // listing all non-special frames
    if(!strstr(frame[i].name,"reload")){
      frame[i].fq/=maxstp;
      frame[i].fq*=N;
      if(frame[i].delay>0){
	frame[i].delay/=maxstp;
	frame[i].delay*=N;
      }
      if(frame[i].idle>0){
	frame[i].idle/=maxstp;
	frame[i].idle*=N;
      }
    }
  }*/

  eprintf("New step / frequence multiplier: %f / (%ld/%ld)\n",dt,N,maxstp);  
  
  delete arr;
  return dt;
}


int PlasmaRecord::Query(){
  int r=Record::Query();
  if(r<0)return r;

  return frame[r].spec;
}


long PlasmaRecord::ReloadGas(Plasma &Gas, int fileadj){
  // test reload frames here

  if(!(wtype&COORD)){
    msg_error("PlasmaRecord.ReloadGas: "
	      "the file '%s' contains no Coord data!\n",file);
    return -1;
  }
  if(!(wtype&VEL)){
    msg_error("PlasmaRecord.ReloadGas: "
	      "the file '%s' contains no Velocity data!\n",file);
    return -1;
  }

  long stp=getNSteps()-1;

  // finding coord and vel frames
  int ivel, icoord, evel, ecoord;
  ivel=FrameWith("ion_veloc");
  icoord=FrameWith("ion_coord");
  evel=FrameWith("electron_veloc");
  ecoord=FrameWith("electron_coord");
  
  if(ivel<0 || icoord<0 || evel <0 || ecoord<0){
    msg_error("PlasmaRecord.ReloadGas: "
	      "the file '%s' contains no sufficient reload data!\n",file);
    return -1;
  }

  
  for(;stp>=0;stp--){
    data_state(stp,1);
    if(frame[ivel].state[1]!=FR_WRITE)continue; // no such frame on this step
    if(frame[icoord].state[1]!=FR_WRITE)continue;
    if(frame[evel].state[1]!=FR_WRITE)continue;
    if(frame[ecoord].state[1]!=FR_WRITE)continue;
    // found !
    eprintf("Trajectory is beeing reloaded from step %ld of %ld.\n",
	    stp,nsteps);

    // reading Gas
    int res=0;
    res+=Get(stp,ivel,Gas.v);
    res+=Get(stp,evel,Gas.v+Gas.ni);
    res+=Get(stp,icoord,Gas.x);
    res+=Get(stp,ecoord,Gas.x+Gas.ni);

    if(res!=4){
      msg_error("PlasmaRecord.ReloadGas: "
		"failed to read particle data from file '%s'\n",file);
      return -1;
    }

    // adjusting file
    if(fileadj){
      Truncate(stp+1);
      step=stp; 
      
    }
    real_step=(step+1e-5)/dstep;  
    
    return step;
  }

  msg_error("PlasmaRecord.ReloadGas: "
	    "cannot find step to reload from in file '%s'\n",file);
  return -1;
}


extern int non_symm;

int PlasmaRecord::AllocPlasma(){
  if(plasmap)delete plasmap;
  plasmap= new Plasma;
  if(!plasmap){
    msg_error("PlasmaRecord: Error allocating Gas!\n");
    return 0;
  }
  plasmap->init(wni,wq,wmass);
  non_symm=plasmap->non_symm=wnsymm;
  plasmap->adjustTG(wT,wGamma);
  return 1;
}



int PlasmaRecord::GetAllNext(){
  long pos=-1;
  int status=0;
  int i;

  for(i=0;i<nfr;i++){
    if(frame[i].state[1]==FR_WRITE){
      if(pos<0)pos=step_position(lastst[1]);
	
      switch(frame[i].spec){
      case PR_IONVEL:
	wake_FP();
	fseek(recfp,pos,SEEK_SET);
	if(fread(plasmap->v,1,frame[i].size,recfp)!=frame[i].size)
	  return status; // file read failure
	status|=frame[i].spec;
	break;
      case PR_ELCVEL:
	wake_FP();
	fseek(recfp,pos,SEEK_SET);
	if(fread(plasmap->v+plasmap->ni,1,frame[i].size,recfp)!=frame[i].size)
	  return status; // file read failure
	status|=frame[i].spec;
	break;
      case PR_IONCOORDS:
	wake_FP();
	fseek(recfp,pos,SEEK_SET);
	if(fread(plasmap->x,1,frame[i].size,recfp)!=frame[i].size)
	  return status; // file read failure
	status|=frame[i].spec;
	break;
      case PR_ELCCOORDS:
	wake_FP();
	fseek(recfp,pos,SEEK_SET);
	if(fread(plasmap->x+plasmap->ni,1,frame[i].size,recfp)!=frame[i].size)
	  return status; // file read failure
	status|=frame[i].spec;
	break;	
      }
	
      pos+=frame[i].size;
      
    }
  }
  return status;
}



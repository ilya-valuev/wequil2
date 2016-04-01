# include "plasma.h"
# include "interact.h"
# include "statist.h"



// units, constants and coefficients

double unit_l;
double RDebye;





Distribution DRRee,DRRep,DRRpp;

/*
double gaussrand(double length,double delta){
  double tmp,val,x;
  do{
    x=2*length*random1();
    val=exp(-x*x/delta);
    tmp=0.5+random1();
    //printf("%f  %f\n",tmp,val);
  }while(tmp>val);
  return x;
}


void maxwellrand(Vector_3 &v,double v_sq){
  double xi[4];
  int i;
  for(i=0;i<4;i++)xi[i]=0.5+random1();

  double a1=sqrt(fabs(-2.*v_sq*log(fabs(xi[0]))));

  v[0]=a1*cos(2*M_PI*xi[2]);
  v[1]=a1*sin(2*M_PI*xi[2]);
  v[2]=sqrt(fabs(-2.*v_sq*log(fabs(xi[1]))))*cos(2*M_PI*xi[3]);
}
*/





void Plasma::fileout(FILE *fp){
  fprintf(fp,"Ne: %d\n",ne);
  fprintf(fp,"Ni: %d\n",ni);
  fprintf(fp,"Ion mass: %f\n",mass);
  fprintf(fp,"Ion charge: %f\n",q);
  //fprintf(fp,"Symmetric: %s\n");
  fprintf(fp,"T: %f (K)\n",par_T*1e4);
  fprintf(fp,"Gamma: %f\n",par_Gamma);
  fprintf(fp,"Electron density: %e (cm^-3)\n",par_density*1e19);
  fprintf(fp,"Wpe: %f (i.u.)\n",Wpe);
  fprintf(fp,"Potential: %s\n",charpot);
  if(non_symm)fprintf(fp,"Symmetric: no\n");
  else fprintf(fp,"Symmetric: yes\n");
}
  



void Plasma::adjustTG(double T,double Gamma){
  par_T=T;
  par_Gamma=Gamma;
  L=pow(4.*M_PI*ne/3,1./3.)/Gamma;
  Wpe=sqrt(4.*M_PI*ne/(L*L*L));
  r0=0.05;
  unit_l=k00/(T*1e4);
  par_density=ne/pow(L*unit_l,3)*1e-25;
  RDebye=1./sqrt(4*M_PI*ne/(L*L*L));
}

void Plasma::adjustGn(double Gamma,double dens){
  par_Gamma=Gamma;
  par_density=dens;
  par_T=pow(4*M_PI*dens*1e25/3, 1./3.)*k00/Gamma*1e-4;
  L=pow(4.*M_PI*ne/3,1./3.)/Gamma;
  Wpe=sqrt(4.*M_PI*ne/(L*L*L));
  r0=0.05;
  unit_l=k00/(par_T*1e4);
  RDebye=1./sqrt(4*M_PI*ne/(L*L*L));
}

void Plasma::adjustTn(double T,double dens){
  par_T=T;
  par_density=dens;
  par_Gamma=pow(4*M_PI*dens*1e25/3, 1./3.)*k00/(T*1e4);
  L=pow(4.*M_PI*ne/3,1./3.)/par_Gamma;
  Wpe=sqrt(4.*M_PI*ne/(L*L*L));
  r0=0.05;
  unit_l=k00/(T*1e4);
  RDebye=1./sqrt(4*M_PI*ne/(L*L*L));
}


// scales velocities of electrons and ions simultaneously
// in symmetric case if T0==T0e considers particles as one big set
// in stable_ions case sets the temperature of electrons to T0, T0e is not used
// if T0e<0 or omitted sets T0e=T0
// T0 is ion temperature

void Plasma::vel_scale(double T0, double T0e){
 if(T0e<0)T0e=T0;
 int i;
 Vector_3 v0e(0,0,0), v0i(0,0,0);

 if(!stable_ions){
  for(i=0;i<ni;i++)v0i+=v[i];
  v0i*=1./ni;

  for(i=ni;i<n;i++)v0e+=v[i];
  v0e*=1./ne;
 }


 Vector_3 vm=(ne*v0e+ni*mass*v0i)/(ne+ni*mass);

 double vsce=0, vsci=0;

 if(!stable_ions){
  for(i=0;i<n;i++){
   v[i]-=vm;
   if(i<ni)vsci+=v[i]*v[i]*mass;
   else vsce+=v[i]*v[i];
  }
 }
 else{
  for(i=ni;i<n;i++){
   vsce+=v[i]*v[i];
  }
 }

 if(!stable_ions){

  if(!non_symm && T0==T0e){
   double k=vsci/vsce;
   Te=2*T0/(k+1);
   Ti=k*Te;

   vsce=sqrt(T0*3.*n/(vsce+vsci));
   vsci=vsce;
  }
  else{
   vsce=sqrt(T0e*3.*ne/vsce);
   vsci=sqrt(T0*3.*ni/vsci);

   Te=T0e;
   Ti=T0;
  }

  for(i=0;i<ni;i++)v[i]*=vsci;
  for(i=ni;i<n;i++)v[i]*=vsce;
 }
 else{
  vsce=sqrt(T0*3.*ne/vsce);
  for(i=ni;i<n;i++)v[i]*=vsce;
  Te=T0;
  Ti=0.;
 }

 if(stable_ions)T=T0;
 else T=(T0*ni+T0e*ne)/(ne+ni);
}


// changes temperature of ions
// does not affect electron distribution

void Plasma::ivel_scale(double T0){
 if(stable_ions){
  Ti=0;
  return;
 }

 int i;
 Vector_3  v0i(0,0,0),v0e(0,0,0);


 //for(i=ni;i<n;i++)v0e+=v[i];
 //v0e*=1./ne;
 for(i=0;i<ni;i++)v0i+=v[i];
 v0i*=1./ni;



 //Vector_3 vm=v0e;
 Vector_3 v_eff;

 double  vsci=0;

 for(i=0;i<ni;i++){
  v[i]-=v0i;
  //v_eff=v[i]-vm;
  v_eff=v[i];
  vsci+=v_eff*v_eff*mass;
 }


 vsci=sqrt(T0*3.*ni/vsci);
 //for(i=0;i<ni;i++)v[i]=vm+(v[i]-vm)*vsci;
 for(i=0;i<ni;i++)v[i]*=vsci;

 Ti=T0;
}

void Plasma::evel_scale(double T0){
 int i;
 Vector_3 v0e(0,0,0), v0i(0,0,0);

 if(!stable_ions){
  for(i=0;i<ni;i++)v0i+=v[i];
  v0i*=1./ni;
  for(i=ni;i<n;i++)v0e+=v[i];
  v0e*=1./ne;
 }


 Vector_3 vm(0,0,0), v_eff;
 if(!stable_ions)vm=v0i;

 double vsce=0;

 for(i=ni;i<n;i++){
  v[i]-=v0e;
  v_eff=v[i]-vm;
  vsce+=v_eff*v_eff;
 }

 vsce=sqrt(T0*3.*ne/vsce);
 for(i=ni;i<n;i++)v[i]=vm+(v[i]-vm)*vsce;

 Te=T0;
}




int Plasma::init_vel(int ndistr=ZEROVEL){
  int i,j;
  int distr;
  double delta,T0;
  for(i=0;i<n;i++){

   if(ndistr==SEPARATE){
    if(i<ni){
     distr=idistr;
     T0=ini_Ti/(rel_Ti ? 1: par_T);
    }
    else{
     distr=edistr;
     T0=ini_Te/(rel_Te ? 1: par_T);
    }
   }
   else{
     distr=ndistr;
     T0=1.;
   }

   for(j=0;j<3;j++){
    if((i<ni && stable_ions) || distr==ZEROVEL)v[i][j]=0;
    else if(distr==MAXWELL_P){
     if(i<ni)delta=T0/mass;
     else delta=T0;
     maxwellrand(v[i].v,delta);
    }
    else if(distr==MAXWELL){
     if(i<ni)delta=2.*T0/mass;
     else delta=2.*T0;
     v[i][j]=gaussrand(10*sqrt(delta),delta);
    }
    else{
      msg_error("Plasma.init_vel : Unknown distribution.\n");
      return 0;
    }
   }
 }

 if(ndistr==SEPARATE){
  if(idistr==MAXWELL && edistr==MAXWELL){
   vel_scale(ini_Ti/(rel_Ti ? 1: par_T),ini_Te/(rel_Te ? 1: par_T)); // simultaneous scale ( better )
  }
  else{
   if(idistr==MAXWELL)ivel_scale(ini_Ti/(rel_Ti ? 1: par_T));
   if(edistr==MAXWELL)evel_scale(ini_Te/(rel_Te ? 1: par_T));
  }
 }
 else if(ndistr==MAXWELL)vel_scale(1.);

 return 1;
}



int Plasma::init_config(double in_cs=-1.,int distr=ZEROVEL){
 int i,j,k;

 if(in_cs<0)in_cs=1.;

 Vector_3 r;
 for(i=0;i<n;i++){// checking distances
   int fertig;
   do{
     for(j=0;j<3;j++){
       x[i][j]=in_cs*L*random1(); //2*L*(0.5-(double)random()/RAND_MAX);

       //if((i<ni && stable_ions) || distr!=MAXWELL)v[i][j]=0;
       //else{
       // double delta;
       // if(i<ni)delta=2./mass;
       // else delta=2.;
       //v[i][j]=gaussrand(sqrt(10.*delta),delta);
       //}
     }
     fertig=1;
     for(j=0;j<i;j++){
       r=x[i]-x[j];
       double R=0;
       for(k=0;k<3;k++){
        if(r[k]>0.5*L)R+=(L-r[k])*(L-r[k]);
        else R+=r[k]*r[k];
       }
       if(R<r0*r0){
         printf("Got: %f\n",sqrt(R));
         fertig=0;
	        break;
       }
     }
   }while(!fertig);
 }

 init_vel(distr);
 //if(distr==MAXWELL)vel_scale(1.);
 return 0;
}


Vector_3  vav(0,0,0),fav(0,0,0), vavi(0,0,0), vave(0,0,0);

int Plasma::stepmove(){
 register int i;
 Vector_3 force;
 v_2=0;
 v_2i=0.;
 vavi=vav=fav=Vector_3(0,0,0);
 if(!stable_ions){
  for(i=0;i<ni;i++){  // ions move
   force=f[i]+ext_force(v[i]);
   v[i]+=force*dt/mass;
   x[i]+=v[i]*dt;
   v_2+=(v[i]*v[i])*mass;
   vav+=v[i]*mass;
   fav+=force;
  }
  v_2i=v_2;
  vavi=vav;
 }
 for(i=ni;i<n;i++){ // electrons move
  force=f[i]+ext_force(v[i]);
  v[i]+=force*dt;
  x[i]+=v[i]*dt;
  v_2+=v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2];
  vav+=v[i];
  fav+=force;
 }
 v_2e=v_2-v_2i;
 vave=vav-vavi;

 if(!stable_ions){
  vav/=(ne+ni*mass); // center-of-mass momentum
  if(one_center){
    Te=(v_2e/ne-2*vav*vave/ne+vav*vav)/3.; // electron temperature
    Ti=(v_2i/ni-2*vav*vavi/ni+vav*vav)/3.; // ion temperature
    T=(v_2/n-vav*vav*(ne+ni*mass)/n)/3.; //temperature
  }
  else{
    Te=(v_2e-vave*vave/ne)/(3*ne);
    Ti=(v_2i-vavi*vavi/(ni*mass))/(3*ni);
    T=(ne*Te+ni*Ti)/n;
  }
 }
 else{
  vav/=ne;
  T=(v_2/ne)/3.;   /// ATTENTION TO THIS !!! -- no average electron velocity counted!!!
  //Etot=1.5*T*ne+Epotent;

  Te=T;
  Ti=0;
 }
 Ekin=0.5*v_2;
 Etot=Ekin /*1.5*T*n*/+Epotent;

 return 0;
}

/*
f0,v0,x0
dtmax,dtmin

amax,damax


int Plasma::propagate(double delta_t){
  double dtleft=delta_t;
  double tmpdt=dt, dtsync;
  lList<int> list0, list1;
  do{
    dt=min(dtmax,dtleft);
    dtsync=dt;

    for(i=0;i<n;i++){
      f0[i]=f[i];
      v0[i]=v[i];
      x0[i]=x[i];
    }
    interaction();
    stepmove();

    inspect_forces(&list0,&list1);
    dtcur=0;
    do{
      while(list1.Count()>0){
        double dta=dt/2;
        if(dta<dtmin){
          eprintf("Minimum step reached with unsatisfactory forces (1)!\n");
          break;
        }
        step_back(dta,list0,list1);
        dt=dta;
        inspect_forces(&list0,&list1);
      }
      dtcur+=dt;
      syncronize(dt,dtsync,list0,list1); //moving initial state to dt, making step forwards with list1
      inspect_forces();
      if(list1.Count>0 && dtsync-dtcur<dtmin){
        eprintf("Minimum step reached with unsatisfactory forces (2)!\n");
        break;
      }

    }while(list1.Count()>0);

    dtleft-=dsync;
  }
  dt=tmpdt;
  return 1;
}


int Plasma::inspect_forces(lList<int> &list0,lList<int> &list1){
  list0.Remove_all();
  list1.Remove_all();
  int i;
  double imass=mass;
  for(i=0;i<n;i++){
    if(i==ni)imass=1;
    if(fabs(v[i]-v0[i])>dvmax || fabs((f[i]-f0[i])/imass)>damax){
      list1.Append(i);
    }
    else list0.Append(i);
  }
  return 1;
}


for(i=0;i<n;i++){
    x[i]=x0[i]+dx[i]*dti;
    v[i]=v0[i]+dv[i]*dti;
    f[i]=f0[i]+df[i]*dti


int step_back(double dti,lList<int> &list0,lList<int> &list1){
  double sh=dti/dt;
  int i;
  Vector_3 tmpf;
  // restoring the state
  for(i=0;i<n;i++){
    x[i]=x0[i];
    v[i]=v0[i];
    f[i]=f0[i];
  }
  double dmpdt=dt;
  stepmove(dti);

  // interpolating forces of list0
  list0.Rewind();
  while(!list0.Done()){
    i=list0.GetCur();
    list0.Next();
    f[i]=f0[i]+df[i]*dti;
  }

  // recalculating forces of list1
  interaction_list();

}

int interaction_list(){
  int i, j;
   // restoring states of list1
  list1.Rewind();
  while(!list1.Done()){
    if(list){
      if(list->Done())break;
      i=list->GetCur();
      list->Next();
    }
    else i++;
    x[i]=x0[i];
    v[i]=v0[i];
    f[i]=f0[i];
  }
  // calculating forces of list
*/




Vector_3 Plasma::getCMVel(int comp){
  int i;
  Vector_3 vi(0,0,0), ve(0,0,0);

  if(comp<0 || comp==0){  // ions
    if(!stable_ions){
      for(i=0;i<ni;i++){
        vi+=v[i];
      }
    }
  }
  if(comp==0)return vi/ni;
  if(comp<0 || comp==1){ // electrons
    for(i=ni;i<n;i++){
      ve+=v[i];
    }
  }
  if(comp==1)return ve/ne;
  // full center-of-mass
  return (vi*mass+ve)/(ni*mass+ne);
}


double Plasma::getT(){
 int i;
 v_2i=v_2=0;
 vavi=vav=Vector_3(0,0,0);

 if(!stable_ions){
  for(i=0;i<ni;i++){
   v_2+=v[i]*v[i]*mass;
   vav+=v[i]*mass;
  }
  vavi=vav;
  v_2i=v_2;
 }
 for(i=ni;i<n;i++){
  v_2+=v[i]*v[i];
  vav+=v[i];
 }
 vave=vav-vavi;
 v_2e=v_2-v_2i;

 if(!stable_ions){
  vav/=(ne+ni*mass); // center-of-mass velocity
  if(one_center){
    Te=(v_2e/ne-2*vav*vave/ne+vav*vav)/3.; // electron temperature
    Ti=(v_2i/ni-2*vav*vavi/ni+vav*vav)/3.; // ion temperature
    T=(v_2/n-vav*vav*(ne+ni*mass)/n)/3.; //temperature
  }
  else{
    Te=(v_2e-vave*vave/ne)/(3*ne);
    Ti=(v_2i-vavi*vavi/(ni*mass))/(3*ni);
    T=(ne*Te+ni*Ti)/n;
  }
 }
 else{
  vav/=ne;
  T=(v_2/ne)/3.;   /// ATTENTION TO THIS !!! -- no average electron velocity counted!!!
  //Etot=1.5*T*ne+Epotent;

  Te=T;
  Ti=0;
 }
 Ekin=0.5*v_2;
 Etot=Ekin /*1.5*T*n*/+Epotent;

 return T;
}


double amod(double a, double b){
  if(a>0){
   while(a>b)a-=b;
  }
  else
    while(a<-b)a+=b;
  return a;
}

double amod1(double a, double b){
  return a-((long)(a/b)*b);
}




Vector_3 null_vect(0,0,0);

# pragma argsused
Vector_3 void_force1(Vector_3& v){ return null_vect;};

# pragma argsused
double void_force2(int type, double R,double &e1, double &e2){
   return 0.;
};




Vector_3 Plasma::rcell(Vector_3 r){
  Vector_3 v;
  int k;
  for(k=0;k<3;k++){
      if(r[k]>0)v[k]=amod1(r[k]+L/2,L)-L/2;
      else v[k]=amod1(r[k]-L/2,L)+L/2;
  }
  return v;
}



int write_distr=0;
char distrfile[256];

//double fucktmp[500];

void Plasma::rcell(int m=-1){
  int i,k;
  for(i=(m<0 ? 0 : m);i<n;i++){
    for(k=0;k<3;k++){
      f[i][k]=0;   // reducing to elementary cell

      if(x[i][k]>0)xx[i][k]=amod1(x[i][k]+L/2,L)-L/2;
      else {
        /*if(x[i][k]<-100000){
          i++;
        } */
        xx[i][k]=amod1(x[i][k]-L/2,L)+L/2;
      }

    }
    if(m>=0)break;
  }
}



int Plasma::interaction(int m, int sum){
 register int i, j,k;
 static Vector_3 r;
 //static Vector_3 df;
 int type=0; // ion-ion

 for(i=(m<0 ? 0 : m);i<n;i++){
   for(k=0;k<3;k++){
     f[i][k]=0;   // reducing to elementary cell

     if(x[i][k]>0)xx[i][k]=amod1(x[i][k]+L/2,L)-L/2;
     else xx[i][k]=amod1(x[i][k]-L/2,L)+L/2;

   }
   if(m>=0)break;
 }


 Quant=Ecoul=0;
 double dEcoul,dEpotent,dQuant;
 double df;

 if(!is_matr || m<0)Epotent=0;


 for(i=0;i<(m<0 ? n : m+1);i++){

   if(i>=ni)type|=0x1;// setting electron-? interaction type
   type&=0x1;         // setting ?-ion interaction type

   for(j=((m<0 || i==m) ? i+1: m); j<n;j++){
    if(j>=ni)type|=0x2;// setting ?-electron interaction type
    //r=xx[i]-xx[j];

    for(k=0;k<3;k++){ // determining the closest
                      //distance and correspondent direction

      r[k]=xx[i][k]-xx[j][k];
      if(r[k]>L/2)r[k]-=L;
      if(r[k]<-L/2)r[k]+=L;
    }

    double R=r.norm();

    if(R<1e-20){
      printf("Got small distance (%d,%d) !\n",i,j);
    }
    r/=R;

    dEcoul=-1/R;
    if(type==ELCELC || type ==IONION){
      dEcoul=-dEcoul;
      if(write_distr){
	if(non_symm && type==IONION)DRRpp.point(R,1.);
	else DRRee.point(R,1.);
      }
    }
    else if(write_distr)DRRep.point(R,1.);

    df=potential(type,R,dEpotent,dQuant);
    ///df=PotentialKELBG(type,R,dEpotent,dQuant);


    //f[i]+=df;
    //f[j]-=df;

    for(k=0;k<3;k++){ // to avoid vector copying
      f[i][k]+=df*r[k];
      f[j][k]-=df*r[k];
    }

    Ecoul+=dEcoul;

    Quant+=dQuant;


    if(is_matr && m>=0)Epotent+=dEpotent-(*umatr)(i,j);
    else Epotent+=dEpotent;

    if(is_matr){
      if(m>=0){
	int l=(i< m ? i : j);
	(*umatr)(l,l)=(*umatr)(i,j); //(l,l) used for storing temporary
	//fucktmp[l]=(*umatr)(i,j);
	//printf("%d <- (%d, %d)[%f]\n",l,i,j,fucktmp[l]);
      }
      (*umatr)(i,j)=dEpotent;
    }

    if(m>=0 && i!=m)break;

   }

 }

 if(is_matr && m>=0 && sum){
   Epotent=0;
   for(i=0;i<n;i++){
     for(j=i+1;j<n;j++){
       Epotent+=(*umatr)(i,j);
     }
   }
 }

 return 0;
}



int PlasmaBunch::interaction(int m, int sum){
 register int i, j,k;
 static Vector_3 r;

 int ret=Plasma::interaction(m,sum);

 // now interacting with bunch particles
 int type=0; // ion-ion
 double dEcoul,dEpotent,dQuant;
 double df;

 for(i=0;i<(m<0 ? n : m+1);i++){

   if(i>=ni)type|=0x1;// setting electron-? interaction type
   else type&=0x2;// setting ion-? interaction type
   if(is_ion)type&=0x1;         // setting ?-ion interaction type
   else  type|=0x2;// setting ?-electron interaction type

   for(j=0;j<nb;j++){

     for(k=0;k<3;k++){ // determining the closest
                       //distance and correspondent direction
       r[k]=xx[i][k]-xb[j][k];
       if(r[k]>L/2)r[k]-=L;
       if(r[k]<-L/2)r[k]+=L;
     }

     double R=r.norm();
     if(R<1e-20)printf("Got small distance to bunch (%d,b%d) !\n",i,j);
     r/=R;

     dEcoul=qb/R;

     df=potential(type,R,dEpotent,dQuant);

     for(k=0;k<3;k++){ // to avoid vector copying
       f[i][k]+=df*r[k];
       //f[j][k]-=df*r[k];
     }
     Ecoul+=dEcoul;
     Quant+=dQuant;
     Epotent+=dEpotent;
   }

 }
 return ret;
}


void PlasmaBunch::propagate_bunch(){
  int i;
  double actdn=dnb+dnrest;
  for(i=0;i<nb;i++){
    xb[i][0]+=vb*dt;
    if(xb[i][0]>=L/2){ // goes out of the cell
      create_particle(i);
      actdn-=1;
    }
  }
  if(actdn<=1)dnrest=actdn;
  else{
    int ncr=(int)actdn;
    dnrest=actdn-ncr;
    if(nb+ncr>limnb){
      msg_error("Too many bunch particles: %d, adding %d?\n",nb,ncr);
      return;
    }
    for(i=0;i<ncr;i++){
      create_particle(nb+i);
    }
    nb+=ncr;
  }
}

void PlasmaBunch::start_bunch(){
  dnb=dt*vb*ratio*ne/L; // how many particles must be born each step
  dnrest=0;
}





# ifndef PLASMA_H
# define PLASMA_H

# include <stdio.h>
# include "common.h"
# include "vector_3.h"
# include "statist.h"

typedef Vector_3 (*fncextfrc)(Vector_3 &);

typedef double   (*fncpotent)(int,double,double &,double&);


extern Vector_3 null_vect;
Vector_3 void_force1(Vector_3& v);
double void_force2(int type, double R,double &e1, double &e2);



extern Distribution DRRee,DRRep,DRRpp;
extern int write_distr;
extern char distrfile[256];

extern double unit_l;
extern double RDebye;


enum dtypes {
 NOTHING =0,
 COORD   =0x1,
 VEL     =0x2,
 FLOW    =0x4 };




# define SEPARATE  3
# define MAXWELL_P 2
# define MAXWELL   1
# define ZEROVEL   0




class Plasma{
 void clear_data(){
   if(x)delete[] x;
   if(xx)delete[] xx;
   if(v)delete[] v;
   if(f)delete[] f;

   if(x0)delete[] x0;
   if(v0)delete[] v0;
   if(f0)delete[] f0;

   if(is_matr)delete umatr;
 }
 int alloc;
public:
 int ne;
 int ni;
 int n;
 Vector_3 *x;
 Vector_3 *xx;// coords reduced to elementary cell
 Vector_3 *v;
 Vector_3 *f;


 Vector_3 *x0, *v0, *f0;  // state storage for variable step scheme
 double mass;
 double q;
 double dt;// time step
 double L; // cell length
 double Ecoul; //Coulomb energy
 double Epotent; // real potential energy
 double Quant; // estimate energy of the corresponding quantum ensemble
 double Etot; // total energy
 double Wpe; // electron plasma frequency

 double v_2;
 double T; // temperature
 double Ekin;  // full kinetic energy, including center-of-mass energy

 double v_2e;
 double Te; // electron temperature
 double v_2i;
 double Ti; // ion temperature
 double ini_Te;  // initial temperatures, distributions,  used by setup
 double ini_Ti;
 int rel_Ti; //  if 1 , the temperatures are relative to par_T
 int rel_Te; //

 int idistr;
 int edistr;


 double r0; // the closest acceptable initial distance


 double par_T;  // temperature set
 double par_Gamma; // Gamma set
 double par_density; // density set
 char dataname[50];

 SimmMatr *umatr;
 int is_matr;

 int stable_ions;
 int non_symm;
 int one_center;   // if 1, uses one common center of mass for temperature determination

 // e^2*Z/(4*pi*eps0)
 double k00;

 Vector_3 (*ext_force)(Vector_3& v); //defines external force
                            // acting on the single particle
 double (*potential)(int type,double R,double &dEp,double &dQuant);
  // defines pair potential
  //and corresponding energy contributions


 char charpot[50];

 Plasma(int nion,int ion_index, double imass=1.){
   if(!init(nion,ion_index,imass)){
     serror("Plasma: Cannot create Plasma object!\n");
   }
 }

 int init(int nion,int ion_index, double imass=1.){
  ni=nion;
  ne=ni*ion_index;
  n=ni+ne;
  x= new Vector_3[n];
  xx= new Vector_3[n];
  v= new Vector_3[n];
  f= new Vector_3[n];

  x0= new Vector_3[n];
  v0= new Vector_3[n];
  f0= new Vector_3[n];


  if(!x || !xx || !f ||!v || !x0 || !f0 ||!v0){
    msg_error("Error allocating memory for data arrays\n");
    return 0;
  }
  alloc=1;
  mass=imass;
  q=ion_index;
  ext_force=(fncextfrc)void_force1;
  potential=(fncpotent)void_force2;
  Ekin=T=Epotent=Ecoul=Quant=Etot=0.;
  is_matr=0;
  umatr=NULL;
  stable_ions=0;
  non_symm=1;
  Te=Ti=0.;
  ini_Te=ini_Ti=1.;
  rel_Te=rel_Ti=1;
  idistr=edistr=0;
  one_center=1;

  //     e^2*Z/k
  k00=(1.602e-19)*(1.602e-19)/(4*M_PI*1.38e-23*8.854e-12);

  return 1;
 }

 Plasma(){  // ONly for storage !!!!
  x=xx=v=f=NULL;
  n=ni=ne=0;
  alloc=0;
  is_matr=0;
  umatr=NULL;
  stable_ions=0;
 }

 ~Plasma(){
   if(alloc){
     clear_data();
   }

 }

 void init_matrix(){
   if(umatr)delete umatr;
   umatr=new SimmMatr(n);
   is_matr=1;
 }


 Plasma &operator=(Plasma &c){
   if(this==&c)return c;
   if(this->alloc){
     clear_data();
   }
   memcpy(this,&c,sizeof(Plasma));
   c.alloc=0;
   return *this;
 }


 virtual int stepmove();
 void rcell(int m);
 Vector_3 rcell(Vector_3 r);
 virtual int interaction(int m=-1, int sum=0);
 int init_config(double,int);
 int init_vel(int distr);
 void vel_scale(double T0, double T0e=-1);
 void ivel_scale(double T0);
 void evel_scale(double T0);
 void adjustTG(double T, double Gamma);
 void adjustTn(double T, double dens);
 void adjustGn(double Gamma, double dens);
 double getT();
 Vector_3 getCMVel(int comp=-1); // gets center-of-mass velocity for the specified component
 void fileout(FILE *fp);
};




// plasma with a bunch propagating through it
// in X direction from left to right
// ratio=n_bunch/ne
// vb= velocity of the bunch
// qb  is the charge of bunch particles (measured in +qe)
class PlasmaBunch: public Plasma{
public:
  int qb;
  double vb;
  int limnb; // limiting number of bunch particles
  int nb;    // actual number of bunch particles
  double dnb; // particles born each step
  double dnrest;
  double ratio;
  Vector_3 *xb; //, *vb;
  int is_ion;
  PlasmaBunch(double sratio, double svb, int sqb,
    int nion,int ion_index, double imass=1.):
    Plasma(nion,ion_index,imass){
    qb=sqb;
    ratio=sratio;

    if(qb>0)is_ion=1; // bunch particle type
    else is_ion=0;

    xb=NULL;
    vb=svb;
    // determining approximate number of bunch particles
    nb=0;
    dnb=0;
    dnrest=0;
    limnb=(int)(2*ne*ratio);
    if(!limnb)return;
    xb = new Vector_3[limnb];
    //vb = new Vector_3[limnb];
    if(!xb /*|| !vb*/){
      fatal_error("PlasmaBunch: constructor MAE!\n");
    }
  }
  ~PlasmaBunch(){
    if(xb) delete [] xb;
    //if(vb) delete [] vb;
  }
  virtual int interaction(int m=-1, int sum=0);
  virtual int stepmove(){
    propagate_bunch();
    return Plasma::stepmove();
  }

  // starts the bunch propagation
  void start_bunch();

  void stop_bunch(){
    nb=0;
    dnb=0;
    dnrest=0;
  }


  void propagate_bunch();

  void create_particle(int i){
    xb[i][0]=-L/2;
    xb[i][1]=L*random1();
    xb[i][2]=L*random1();
  }

};


# endif





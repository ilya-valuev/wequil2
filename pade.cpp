# include <stdio.h>
# include <math.h>
# include "common.h"



#define SQR3 1.7320508
#define CE   0.577216
#define LN2  0.69314718
#define LN3  1.0986123

double a1=10.134;
double a2=9.8227;

double u0=-0.5*SQR3;
double u1=-9./8.;
double u2=-9./8.*LN3-3./2.*CE+1.;
double u3= -27./16.*SQR3;
double u4= -9./4.*SQR3*(1./4.*LN3+LN2+CE-5./4.-3./(8*M_PI*M_PI)*(a1-a2/(4*M_PI)));

double U2mc(double G){
  double a=-0.893831;
  double b=0.130345;
  double c=0.321362;
  //double d=0;
  double r=1.00076;
  double s=0.186181;

  double U2_mc=a*pow(G,r)+(b*log(G)+c)*pow(G,s); //+d/pow(G,t);
  return U2_mc;
}

double  Upade(float G){
  double logG=log(G);
  double u03=u0;
  double u13=u1-u0*u3/u2+u0*u1*u4/u2/u2;
  double u23=u2-u0*u4/u2;
  double u33=-u3/u2+u1*u4/u2/u2;
  double u43=-u4/u2;
  double u53=-u1/u0*u33;
  double gamma=1481.83;
  double delta=9.21229;

  double Up=(u03*pow(G,1.5)+(u13*logG+u23)*pow(G,3)+gamma*pow(G,delta)*U2mc(G) )/(1.+(u33*logG+u43)*pow(G,1.5)+u53*pow(G,3.)*logG*logG+gamma*pow(G,delta));

  
  return Up;
}




void main(){
  FILE *f1=fopen("pade.dat","wt");

  double g;

  for(g=0.01;g<50; g+=0.02){
    fprintf(f1,"%e %e\n",g,Upade(g));
  }
  
  fclose(f1);
}

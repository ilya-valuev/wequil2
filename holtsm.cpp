# include <stdio.h>
# include <math.h>
# include <stdlib.h>
# include "common.h"

double Holtsmark(double x){
  if(x<=1e-6)return 0.;

  double dt=0.01*2*M_PI/x;
  double t=dt,sum=0.,ds;
  do{
   ds=exp(-pow(t,3./2.))*sin(x*t)*t*dt;
   sum+=ds;
   t+=dt;
  }while(t<5.);

  return (2/M_PI)*x*sum;
}


void write_holtsmark(){
  FILE *f1=fopen("holtsm.d","wt");
  for(double x=0.; x< 20; x+=0.01){
    fprintf(f1,"%g %g\n",x,Holtsmark(x));
  }
  fclose(f1);
}


double Sum(double xi){
  if(xi<0)return 0.;
  double s=0;
  double i;
  for(i=1;i<100;i++){
    s+=(1.+xi*xi/(4*i*i))/(i*i*i);
  }
  return (double)s*sqrt(M_PI)*xi*xi*xi;
}

double Sum_full(double xi){
  if(xi<0)return 0.;
  double s=0,ds;
  double i;
  for(i=1;i<1000.;i+=1){
    ds=exp(xi*xi/(4*i*i))/(i*i*i);
    s+=ds;
    if(ds<1e-2)break;
  }
  if(i>=1000.)printf("! Accuracy is: %f\n",ds);

  return (double)s*sqrt(M_PI)*xi*xi*xi;
}





double Integ(double xi){
 
  double dt=1e-4;
  double t=dt,sum=0.,ds;
  do{
   ds=t*exp(-t*t)/(1.-exp(-xi*M_PI/t))*dt;
   sum+=ds;
   t+=dt;
   //printf("%e %e %e\n",t,ds,sum);
  }while(t<15 /*fabs(sum)<1e-32 || fabs(ds/sum)>1e-3*/);

  //printf("OK\n");
  return 4*sqrt(M_PI)*xi*sum;
}



main(){
  write_holtsmark();
  exit(0);
  double x=0, dx=0.005, maxx=10;
  

  double me=0.9109534e-30;
  double kb=1.389622e-23;
  double hq=1.0545887e-34;

  double t_k=(1.602e-19)*(1.602e-19)/(4*M_PI*8.854e-12)*sqrt(me)*pow(kb*1e4,-3./2.);
  double h_k=hq/(1e4*kb*t_k);
  double lambda_k=h_k/sqrt(2.);
  double xi_k=1./lambda_k;
  
  printf("t_k: %e\n h_k: %e\nlambda_k: %e\nxi_k: %e\n",t_k,h_k,lambda_k,xi_k);

#if 0
  int i;
  for(i=10;i<=20;i++){
    double S1=Integ(i/M_PI);
    double S2=Sum_full(i/M_PI);
    double S3=0.5*Integ(-i/M_PI);
    printf("%d  %1.4e %1.4e %1.4e %1.4e %1.4e %1.4e\n",i, S1,S2,S1+S2,S3, log(S1+S2),log(S3));
  }
# endif

  double fee, fie;
  FILE *f1=fopen("gnus.gnu","wt");
  for(x=8.;x<maxx;x+=dx){
    double xi=xi_k/sqrt(fabs(x));
    
    
    fee=0.5*Integ(-xi);
    fie=Integ(xi)+Sum(xi);
    //if(f<1e-32)lnf= -32;
    //else if(f>1e32)lnf=32;
    //else lnf=log(f);
    printf("%f\n",x);
    fprintf(f1,"%e %e %e %e %e\n",x,-log(fee),-log(fie),sqrt(M_PI)*xi,-sqrt(M_PI)*xi);
  }
  fclose(f1);
  return 0;
}

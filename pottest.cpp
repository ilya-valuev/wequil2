# include <stdio.h>
# include <math.h>

# include "interact.h"

int main(){

 float T;
 double uee, uei,uq;

 Pauli_part=0;
 equal_lambda=1;
 non_symm=0;

 FILE *f=fopen("test.dat","wt");
 
 for(T=1.;T<=10.;T+=1.){
  Lambda_ee=0.179*sqrt(T);
  Lambda_pp=Lambda_ep=Lambda_pauli=Lambda_ee;

  Correction(IONION,PotentialKELBG,T);
  Correction(IONELC,PotentialKELBG,T);
  Correction(ELCELC,PotentialKELBG,T);

  PotentialKELBG(ELCELC,1e-5,uee,uq);
  PotentialKELBG(IONELC,1e-5,uei,uq);

  fprintf(f,"%f %f %f\n",T*1e4,-uee,-uei);
 }
 fclose(f);
 return 0;
}
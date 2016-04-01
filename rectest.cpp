# include "common.h"

# include <stdio.h>


# include "record.h"


int main(int argc, char *argv[]){
  if(argc<2){
    serror("Usage: rectest <recfile>\n");
  }

  Record Rec;
  Rec.Open(argv[1]);

  int i,n;
  n=Rec.getNFrames();
  char fname[500];
  FILE *f;
  long Nstp=Rec.getNSteps();

  for(i=0;i<n;i++){
    sprintf(fname, "frame%d.dat",i);
    f=Err_fopen(fname,"wt");

    long j,sz;
    
    for(j=0;j<Nstp;j++){
      sz=Rec.FrameState(j,i);
      if(sz){
	fprintf(f,"%ld %ld\n",j,sz);
      }
    }
    fclose(f);
  }
  return 1;
}

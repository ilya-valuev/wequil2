# include <stdio.h>


# include "common.h"
# include "rparam.h"
# include "analyser.h"



int main(int argc, char *argv[]){

# ifndef UNIX
  Exit_wait(1);
# endif

  char *cfgfile="analyse.cfg";
  if(argc>1)cfgfile=argv[1];

  Set_comment_char(';');
  Open_param_file(cfgfile);

  char datafile[256];
  if(argc<=2){
   Read_param("Trajectory file: %s",datafile);
  }
  else strcpy(datafile,argv[2]);

  Close_param_file();

  PlasmaRecord Trj;

  Trj.Open(datafile);

  AnalyseShell AnProc;


  AnProc.Setrec(&Trj);

  AnProc.Init(cfgfile);


  if(argc>=4){
    AnProc.SetOutDir(argv[3]);
  }

  AnProc.StartRead();

  int part=0;
  do{

    AnProc.StepRec();
    if(part<(int)(100.*AnProc.ndone())){
      part++;
      if(part%2==0)printf("%d %% of trajectory processed, record=%d, time=%lf.\n",part,AnProc.cur_rp,AnProc.rtime);
    }

    // printf("hkqwejvfkjq: %d\n",AnProc.trj_end);

  }while(!AnProc.trj_end);

  AnProc.Process();

# ifndef UNIX
  msg_error("Finished!\n");
# endif  
  return 0;
}









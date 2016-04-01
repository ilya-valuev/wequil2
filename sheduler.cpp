# include<stdio.h>


# define MAX_RUN 25

typedef struct {
  char name[20];
  char file[200];
  char host[20];
} execrun;


execrun run[MAX_RUN];

cList actionlist;
NamedList indruns;

char cequil_dir[200], analyse_dir[200];

void cmd_analysecfg(char  *cfgfile){
  actionlist.rewind();
  char name[20];
  int irun;
  do{
    sprintf(name,"%d",actionlist.step());
    irun=indruns.search(name);
    if(irun<0)serror("Run %s -- no such run in the list!\n",name);
    
    Open_param_file(run[irun].file);
    char proc_dir[200], dname[20];
    Read_param("Process output directory: %s",proc_dir);
    Read_param("Data name: %s",dname);
    if(proc_dir[strlen(proc_dir)]!='/')strcat(proc_dir,"/");
    strcat(proc_dir,dname);
    strcat(proc_dir,"1.trj");

    Close_param_file();
    
    char ancfg[200], tmpstr[200];
    strcpy(ancfg,analyse_dir);
    
    sprintf(tmpstr,"analyse.%s.%s",run[irun].host,run[irun].dname);
    strcat(ancfg,tmpstr);

    sprintf(tmpstr,"cp %s %s",cfgfile,ancfg);
    if(system(tmpstr)==-1)printf("\nExec: %s\n",strerror(errno));
    
    Open_param_file(ancfg,"w");
    Write_param("Trajectory file:",proc_dir);
    
    


int main(int argc, char *argv[]){

  Open_param_file("sheduler.cfg");

  Read_param("cequil directory: %s", cequil_dir);
  if(cequil_dir[strlen(cequil_dir)]!='/')strcat(cequil_dir,"/");
  Read_param("analyse directory: %s", analyse_dir);
  if(analyse_dir[strlen(analyse_dir)]!='/')strcat(analyse_dir,"/");

  Set_position("List of runs:");
  int nrun;
  do{
    nrun=0;
    if(Read_param("$:%s %s %s",run[nrun].name,run[nrun].dir,run[nrun].host)!=3){
      if(strstr(AcFormat,"end"))break;
      serror("Invalid run list in sheduler.cfg!\n");
    }
    indruns.insert(AcFormat);
    nrun++;
    if(nrun>=MAX_RUN)printf("Too many runs in the list, skipping from #%d!\n",nrun+1);
  }while(nrun<MAX_RUN];

  
  Set_position("Commands:");
  char cmd_str[200], cmd1[50],cmd2[50],cfgfile[200];
  

  do{
    Read_param("$:>",cmd_str);

    sscanf(cmd_str,"%s",cmd1);
    if(!strcmp(cmdstr,"analysecfg")){
      if(sscanf(cmd_str,"%s %s %s",cmd1,cfgfile,cmd2)<2)
	serror("Command usage: analysecfg cfgfile [-Rrunlist]\n");

      if(!strncmp(cmd2,"-R",2))actionlist=cList(cmd2+2);
      else actionlist=cList(0,nrun);

      Close_param_file();
      
      cmd_analysecfg(cfgfile);
    }
  }
      

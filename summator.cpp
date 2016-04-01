# include <stdio.h>
# include <string.h>

# include "common.h"
# include "rparam.h"
# include "four.h"

int main(int argc,char *argv[]){
  Exit_wait(1);
  int cols, cole, nrep, nlines;
  char dname[250], oname[20];
  char str[500],fname[250];
  char cfgfile[250]="summator.cfg";

  if(argc>1)strcpy(cfgfile,argv[1]);

  Open_param_file(cfgfile);

  if(argc>2)strcpy(dname,argv[2]);
  else Read_param("Data name: %s",dname);

  Read_param("Start column: %d", &cols);
  Read_param("End column: %d", &cole);
  Read_param("Repeats: %d", &nrep);

  if(argc>3)strcpy(oname,argv[3]);
  else Read_param("Output name: %s",oname);

  Read_param("Number of lines: %d",&nlines);

  int olines=nlines;
  if(Read_param("Lines in output: %d",&olines)){
   if(olines>nlines)fatal_error("Number of lines in output is greater than in input!\n");
  }

  Read_param("Fourier: %s",str);
  int make_four=0;
  int tcol=0;
  double Wstart=0,Wend=0,dW=0;
  if(strstr(str,"y"))make_four=1;
  else if(!strcmp(str,"complex"))make_four=2;
  else if(!strcmp(str,"complex_sq"))make_four=3;
  int corr_range=0;
  if(make_four){

   Read_param("Correct range: %s",str);
   if(strstr(str,"y"))corr_range=1;
   else{
    Read_param("Range: %lf, %lf",&Wstart,&Wend);
    dW=(Wend-Wstart)/(olines-1);
   }
  }

  int integ=0;
  Set_stop(0);
  if(Read_param("Integration: %s",str)){
   if(strstr(str,"y"))integ=1;
  }
  int fluc=0;
  char dispfile[256];
  if(Read_param("Dispersion: %s",dispfile)){
   fluc=1;
   if(strstr(dispfile,"default")){
    set_extension(dispfile,oname,"dis");
   }
  }
  int sline=0;
  if(!Read_param("Start line: %d",&sline))sline=0;



  Set_stop(1);
  if(make_four || integ){
   Read_param("Time column: %d",&tcol);
  }
  Close_param_file();
  int i,j,m;

  double *arr;
  int ncols=cole-cols+1;
  arr= new double[2*nlines*ncols];
  if(!arr)serror("MAE\n");

  double *arr2;
  if(fluc){
   arr2= new double[2*nlines*ncols];
   if(!arr2)serror("MAE\n");
  }

  double *ints=NULL;
  if(integ){
   ints= new double[ncols];
   if(!ints)serror("MAE\n");
   for(i=0;i<ncols;i++)ints[i]=0;
  }

  for(i=0;i<2*nlines*ncols;i++)arr[i]=0.;
  if(fluc){
   for(i=0;i<nlines*ncols;i++)arr2[i]=0.;
  }
  double val;
  double time_int=0.;

  FILE *f1;
  for(i=0;i<nrep;i++){

   /*if(i>0)sprintf(str,"%d",i);
   else sprintf(str,"%s","");*/
   sprintf(str,"%d",i);

   sprintf(fname,dname,str);
   f1=Err_fopen(fname,"rt");
   printf("File #%d...\n",i+1);

   for(m=0;m<nlines;m++){


    if(fgetline(f1,str,500)==0){
     if(feof(f1))fatal_error("End of file reached at line #%d!\n",m);
     m--;
     continue; // void line
    }
    if(str[0]=='#'){
     m--;
     continue;
    }


    char *scan=strtok(str," ");
    int k=0;
    j=0;
    while(scan){
      if(i==0 && (make_four || integ)){ // getting time interval
        if(j==tcol){
          if(m==sline)sscanf(scan,"%lf",&time_int); // first t
          else if(m==nlines-1){ // last t
            sscanf(scan,"%lf",&val);
            time_int=fabs(val-time_int);
            if(corr_range){
              Wstart=0;
              dW=2*M_PI/time_int;
              Wend=Wstart+dW*(olines-1);
            }
          }
        }
      }
      if(j>=cols && j<=cole){
       sscanf(scan,"%lf",&val);
       arr[2*nlines*k+2*m]+=val;
       if(fluc){
        arr2[nlines*k+m]+=val*val;
       }
       k++;
      }
      scan=strtok(NULL," ");
      j++;
    }
   }
   fclose(f1);
  }

  for(i=0;i<nlines;i++){
   for(j=0;j<ncols;j++){
     arr[2*nlines*j+2*i]/=nrep;
     if(fluc){ // dispersion
      arr2[nlines*j+i]/=nrep;
      arr2[nlines*j+i]-=arr[2*nlines*j+2*i]*arr[2*nlines*j+2*i];
      arr2[nlines*j+i]=sqrt(arr2[nlines*j+i]);
     }
   }
  }

  if(make_four){
    printf("Fourier transform...\n");
    double sc=time_int/(2*M_PI);

    if(make_four>=2){ // complex transform
     for(j=0;j<ncols;j+=2){
       printf("Columns %d and %d ...\n",cols+j,cols+j+1);
       double t0=time_int/(nlines-1);


       double *tmp;
       for(m=0;m<nlines;m++){
        arr[2*nlines*j+2*m]*=sc;
        arr[2*nlines*j+2*m+1]=sc*arr[2*nlines*(j+1)+2*m]; // copying im part from the next
       }

       tmp=four_direct(arr+2*nlines*j,nlines,1,t0*Wstart,t0*Wend,t0*dW);
       if(make_four==2){
        for(m=0;m<olines;m++){
          arr[2*nlines*j+2*m]=tmp[2*m];
          arr[2*nlines*(j+1)+2*m]=tmp[2*m+1];
         //arr[nlines*(j+1)+2*m]=sqrt(tmp[2*m]*tmp[2*m]+tmp[2*m+1]*tmp[2*m+1]);
        }
       }
       else{
        for(m=0;m<olines;m++){
          arr[2*nlines*j+2*m]=tmp[2*m];
          //arr[2*nlines*(j+1)+2*m]=tmp[2*m+1];
          arr[2*nlines*(j+1)+2*m]=sqrt(tmp[2*m]*tmp[2*m]+tmp[2*m+1]*tmp[2*m+1]);
        }
       }
       delete tmp;
     }
    }
    else{
     for(j=0;j<ncols;j++){
       printf("Column %d ...\n",cols+j);
       double t0=time_int/(nlines-1);
       for(m=0;m<nlines;m++){
        arr[2*nlines*j+2*m]*=sc;
        arr[2*nlines*j+2*m+1]*=sc;
       }

       double *tmp;
       tmp=four_direct(arr+2*nlines*j,nlines,1,t0*Wstart,t0*Wend,t0*dW);
       for(m=0;m<olines;m++){
         arr[2*nlines*j+2*m]=sqrt(tmp[2*m]*tmp[2*m]+tmp[2*m+1]*tmp[2*m+1]);
       }
       delete tmp;
     }
    }
  }



  sprintf(fname,dname,"0");
  f1=Err_fopen(fname,"rt");

  FILE *f2=Err_fopen(oname,"wt");
  FILE *f3;
  if(fluc){
   f3=Err_fopen(dispfile,"wt");
  }
  printf("Writing data...\n");

  for(m=0;m<olines;m++){

   fgetline(f1,str,500);
   if(str[0]=='#'){
    fprintf(f2,"%s\n",str);
    if(fluc)fprintf(f3,"%s\n",str);
    m--;
    continue;
   }

   char *scan=strtok(str," ");
   int k=0;
   j=0;
   while(scan){
     if(make_four && j==tcol){
      fprintf(f2,"%13e ",Wstart+m*dW);
      if(fluc)fprintf(f3,"%13e ",Wstart+m*dW);
     }
     else if(j>=cols && j<=cole){
      fprintf(f2,"%13e ",arr[2*nlines*k+2*m]);
      if(fluc)fprintf(f3,"%13e ",arr2[nlines*k+m]);
      if(integ && m>=sline){
       ints[k]+=arr[2*nlines*k+2*m];
      }
      k++;
     }
     else{
      sscanf(scan,"%lf", &val);
      fprintf(f2,"%13e ",val);
      if(fluc)fprintf(f3,"%13e ",val);
     }
     scan=strtok(NULL," ");
     j++;
   }
   fprintf(f2,"\n");
   if(fluc)fprintf(f3,"\n");
  }
  if(integ){
   fprintf(f2,"#integ: ");
   for(i=0;i<ncols;i++){
    fprintf(f2,"%d:%f ",cols+i+1,ints[i]*time_int/nlines);
   }
   fprintf(f2,"\n");
  }


  fclose(f1);
  fclose(f2);
  if(fluc)fclose(f3);
  delete arr;
  if(fluc)delete arr2;
  printf("Finished.\n");
  //printf("Press ENTER to exit\n");
  //getchar();
  return 0;
}
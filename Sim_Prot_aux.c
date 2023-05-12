#include "Sim_Prot_aux.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

char Get_compression(char *pdb_name){
  char *tmp=pdb_name;
  while((*tmp!='\0')&&(*tmp!='\n')){
    if(*tmp=='.'){
      if((*(tmp+1)=='g')&&(*(tmp+2)=='z')){
	return('g');
      }else if((*(tmp+1)=='Z')&&(*(tmp+2)=='\0')){
	return('Z');
      }else if((*(tmp+1)=='z')&&(*(tmp+2)=='i')&&
	       (*(tmp+3)=='p')&&(*(tmp+4)=='\0')){
	return('z');
      }
    }
    tmp++;
  }
  return('\0');
}

int Num_lines(char *file_name){
  FILE *file_in=fopen(file_name, "r"); int n=0;
  char string[5000];
  if(file_in==NULL){
    printf("WARNING, file %s does not exist\n", file_name);
    return(0);
  }
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]!='#')n++;
  }
  fclose(file_in);
  return(n);
}

int Find_name(char *name, char **names, int N, int Np)
{
  int ip;
  for(ip=Np; ip<N; ip++)if(strcmp(name, names[ip])==0)return(ip);
  for(ip=0; ip<Np; ip++)if(strcmp(name, names[ip])==0)return(ip);
  return(-1);
}

int Change_ext(char *name_out, char *name_in, char *ext){
  char *ptr1=name_in, *ptr2=name_out;
  while((*ptr1!='\0')&&(*ptr1!='.')){
    *ptr2=*ptr1; ptr1++; ptr2++;
  }
  sprintf(name_out, "%s%s", name_in, ext);
  //sprintf(name_out, "%s%s", name_out, ext);
  printf("file name: %s\n", name_out);
  return(0);
}

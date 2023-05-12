#include <stdio.h>
#include <stdlib.h>
#include "Print_pairwise.h"
extern void Change_ext(char *, char *, char *);

void Print_pairwise(int ***Ali_pair, // Ali_pair[i][j][site_i]=site_j i<j
		    //int **PC_opt_pair,
		    int *Lprot,      // Length of protein i
		    int Nprot,       // Number of proteins
		    char **Seq,      // All sequences
		    char **name_seq, // Names of sequences
		    char *name_in)
{
  char name_pwa[100];
  Change_ext(name_pwa, name_in, ".pwa");
  printf("Printing pairwise PC alignments in file %s\n", name_pwa);
  int iprot, jprot, i;
  FILE *file_out=fopen(name_pwa, "w");
  for(iprot=0; iprot<Nprot; iprot++){
    int Li=Lprot[iprot];
    char *Seq_i=Seq[iprot];
    char *name_i=name_seq[iprot];
    for(jprot=0; jprot<iprot; jprot++){
      char *Seq_j=Seq[jprot]; int Lj=Lprot[jprot];
      fprintf(file_out, ">%s_%s %d %d %d %d\n",  //PC_opt= %d
	      name_i, name_seq[jprot], iprot, Li, jprot, Lj);
      //PC_opt_pair[iprot][jprot]
      /*fprintf(file_out, ">%d_%d %d %s L= %d %d %s L= %d\n",  //PC_opt= %d
	      iprot, jprot, iprot, name_i, L,
	      jprot, name_seq[jprot], Lprot[jprot]);*/
      /*for(i=0; i<L; i++)fprintf(file_out, "%c", Seq_i[i]);
      fprintf(file_out, "\n");
      for(i=0; i<L; i++){
	if(ali[i]<0){fprintf(file_out, "-");}
	else{fprintf(file_out, "%c", Seq_j[ali[i]]);}
      }
      fprintf(file_out, "\n"); */
      int *ali=Ali_pair[iprot][jprot];
      //for(i=0; i<Li; i++)fprintf(file_out, "%d ", ali[i]);
      int jj=0; // print sequence i
      for(i=0; i<Li; i++){
	int j=ali[i];
	if(j>=0){
	  while(jj<j){fprintf(file_out, "-"); jj++;} jj++;
	}
	fprintf(file_out, "%c", Seq_i[i]);
      }
      while(jj<Lj){fprintf(file_out, "-"); jj++;}
      fprintf(file_out, "\n");
      jj=0;  // print sequence j
      for(i=0; i<Li; i++){
	int j=ali[i];
	if(j>=0){
	  while(jj<j){fprintf(file_out,"%c", Seq_j[jj]); jj++;} jj++;
	  fprintf(file_out, "%c", Seq_j[j]);
	}else{
	  fprintf(file_out, "-");
	}
      }
      while(jj<Lj){fprintf(file_out, "%c", Seq_j[jj]); jj++;}
      fprintf(file_out, "\n");
    }
  }
  fclose(file_out);
  return;
}

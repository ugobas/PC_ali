#include "Contact_divergence_aux.h"
#include "PC_ali.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define TM_coeff 1.24
#define TM_offset 1.8
#define L_offset 15

void Closest_neighbor_PC(int *neighb1, float *score1,
			 float **score, int n1, int n2, int order);

void Align_PC(int *ali_PC, int *ali, double *PC_load, float d02,
	      float **d2, int **nc, char *seq1, int n1, char *seq2, int n2)
{
  // Compute score
  float *score[n1];
  for(int i=0; i<n1; i++){
    score[i]=malloc(n2*sizeof(float));
    for(int j=0; j<n2; j++){
      double S=0; // no score for aligned pair (it is constant)
      if(seq1[i]==seq2[j])S+=PC_load[1];
      S+=PC_load[2]*d02/(d02+d2[i][j]);
      S+=PC_load[3]*nc[i][j];
      score[i][j]=S;
    }
  }

  float thr=PC_load[1]*0.2+PC_load[2]/(1+0.5)+PC_load[3]*2;
  // Minimum score for alignment
  int neighb1[n1], neighb2[n2]; // ali_tmp[n1];
  float score1[n1], score2[n2];
  Closest_neighbor_PC(neighb1, score1, score, n1, n2, 0);
  Closest_neighbor_PC(neighb2, score2, score, n2, n1, 1);
  Align_neighbors(ali_PC, ali, thr, neighb1, score1, n1, neighb2,score2,n2);
  for(int i=0; i<n1; i++)free(score[i]);
}

void Closest_neighbor_PC(int *neighb1, float *score1,
			 float **score, int n1, int n2, int order)
{
  for(int i=0; i<n1; i++){
    float s_max=0, s; int j_max=-1;
    for(int j=0; j<n2; j++){
      if(order==0){s=score[i][j];}else{s=score[j][i];}
      if(s>s_max){s_max=s; j_max=j;}
    }
    neighb1[i]=j_max; score1[i]=s_max;
  }
}

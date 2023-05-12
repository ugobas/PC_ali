int IT_MAX=2;
float penalty=2.6; // Gap insertion penalty value
float penext=0.1; // Extension penalty
float penalty_ini=3.0; // Gap insertion penalty value
float penext_ini=0.1; // Extension penalty
int NORM_SIM=0;
int GAP_SS=0;

/*float weight[5]= // Optimized on 50044
  {0.75,  // weight of ali
   0.85,  // weight of seq
   0.80,  // weight of sec.str.
   0.90,  // weight of coord (TM)
   0.90   // weight of contacts
   }; */
float weight[5]= 
  {0.3,  // weight of ali
   0.6,  // weight of seq
   0.5,  // weight of sec.str.
   1.0,  // weight of coord (TM)
   1.0   // weight of contacts
   };
int ini_ss=1; // Subtract gap scores from SS scores

// Parameters of the TM score
#define TM_coeff 1.24
#define TM_offset 1.8
#define L_offset 15

/*************************************************************************

   Needleman-Wunsch subroutine based on the Program:    
   align.c
   
   Version:    V3.2
   Date:       27.02.07
   Function:   Perform Needleman & Wunsch sequence alignment
   
   Copyright:  (c) SciTech Software 1993-2007
   Author:     Dr. Andrew C. R. Martin
   EMail:      andrew@bioinf.org.uk
               
**************************************************************************

   This program is not in the public domain, but it may be copied
   according to the conditions laid out in the accompanying file
   COPYING.DOC

   The code may not be sold commercially or included as part of a 
   commercial product except as described in the file COPYING.DOC.

**************************************************************************

   Description:
   ============
   A simple Needleman & Wunsch Dynamic Programming alignment of 2 
   sequences.  
   A window is not used so the routine may be a bit slow on long 
   sequences.

**************************************************************************/

#include "protein.h"
#include "nj_align.h"
#include "allocate.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tm_score.h"
#include "Contact_divergence_aux.h"

//#include "blosum62.h"
//#include "secstrsim.h"
extern int blosum62[22][22];
extern float SS_sim_score[4][4];

void Fill_sim_matrix(float **sim_mat, int **nss1, int **nss2,
		     int **msa, struct protein **prot_p,
		     float **d2, int **nc, int comp_nc,
		     int ***Ali_pair, float ****Rot_pair,
		     float ***Shift_pair,
		     int *index1, int n1, int l1,
		     int *index2, int n2, int l2);
void Sum_ss(int **nss, short *ss3, int *ires, int l);
void Add_aligned(float **sim_mat, float score,
		 int len1, int *ires_1, int l1,
		 int len2, int *ires_2, int l2);
void Add_seq_sim(float **sim_mat, float score,
		 short *seq1, int *ires_1, int l1,
		 short *seq2, int *ires_2, int l2);
void Add_sec_str(float **sim_mat, float score,
		 short *ss1, int *ires_1, int l1,
		 short *ss2, int *ires_2, int l2);
void Add_TM(float **sim_mat, float score, float **d2,
	    float **Rot, float *Shift,
	    struct protein *prot1, int *ires_1, int l1,
	    struct protein *prot2, int *ires_2, int l2);
void Add_CO(float **sim_mat, float score, int **nc, int *ali, int comp_nc,
	    short **Cont_map1, int len1, int *ires_1, int l1,
	    short **Cont_map2, int len2, int *ires_2, int l2);
static void Write_msa(int **msa, int *index, int n, int *ali, int lali);

extern void Rotate_vector(float *x, float **rot);
extern void All_distances(float **d2, float *x1, int n1, float *x2, int n2);

int affinealign_new(float **sim_mat, int **nss1, int **nss2,
		    int  l1, int  l2, 
		    int verbose, float  penalty, float  penext,
		    int *align1, int *align2);

extern int affinealign(float **sim_mat, int  l1, int  l2, 
		       int verbose, float  penalty, float  penext,
		       int *align1, int *align2);

static int  SearchForBest(float **matrix, int len1, int len2, 
                          int *BestI, int *BestJ, int *align1, int *align2);
static float TraceBack(// input:
		       float **matrix,int **dirx,int **diry,int len1,int len2,
		       // output:
		       int *align1, int *align2, int *align_len);

/****************** Main routine ********************************/
int Align_PC_NW_mult(int **msa, int *len_msa,
		     int *N_ali_max, int N_seq,
		     struct protein **prot_p,
		     int ***Ali_pair, float ****Rot_pair, float ***Shift_pair,
		     int *index1, int n1, int *index2, int n2)
{
  int l1=len_msa[index1[0]], l2=len_msa[index2[0]], l1_ini=l1, l2_ini=l2;
  int lali=l1+l2, lmax=lali, vbs=0;

  float **sim_mat=Allocate_mat2_f(*N_ali_max, *N_ali_max);
  int **nss1=NULL, **nss2=NULL;
  if(GAP_SS){
    nss1=Allocate_mat2_i(*N_ali_max, 2);
    nss2=Allocate_mat2_i(*N_ali_max, 2);
    printf("gap(c)=%.4f gap(H)=%.4f gap(E)=%.4f\n",
	   SS_sim_score[0][3], SS_sim_score[1][3], SS_sim_score[2][3]);
  }

  int comp_nc=1;
  float **d2=Allocate_mat2_f(lmax, lmax);
  int   **nc=Allocate_mat2_i(lmax, lmax);
  for(int iter=0; iter<IT_MAX; iter++){
    if(iter==0){ 
      Fill_sim_matrix(sim_mat, nss1, nss2, msa, prot_p, d2, nc, comp_nc,
		      Ali_pair, Rot_pair, Shift_pair,
		      index1, n1, l1, index2, n2, l2);
    }else{
      Fill_sim_matrix(sim_mat, nss1, nss2, msa, prot_p, d2, nc, comp_nc,
		      NULL, NULL, NULL,
		      index1, n1, l1, index2, n2, l2);
    }
    int ali1[l1+l2], ali2[l1+l2]; 
    lali=affinealign_new(sim_mat, nss1, nss2, l1, l2, vbs,
			 penalty, penext, ali1, ali2);

    if(lali > *N_ali_max){
      int Na_new=lali*1.2;
      printf("WARNING in Align_PC, lali= %d > max.ali= %d\n",
	     lali, *N_ali_max);
      printf("Increasing  max.ali to %d\n", Na_new);
      for(int i=0; i<N_seq; i++){
	int tmp[Na_new], *t=tmp, *m=msa[i], j;
	for(j=0; j<Na_new; j++){
	  if(j<*N_ali_max){*t=*m; m++;}else{*t=-1;} t++; 
	}
	free(msa[i]);
	msa[i]=malloc(Na_new*sizeof(int));
	t=tmp; m=msa[i];
	for(j=0; j<Na_new; j++){*m=*t; t++; m++;}
      }
      Empty_matrix_f(sim_mat, *N_ali_max);
      sim_mat=Allocate_mat2_f(Na_new, Na_new);
      if(nss1){
	Empty_matrix_i(nss1, *N_ali_max);
	Empty_matrix_i(nss2, *N_ali_max);
	nss1=Allocate_mat2_i(Na_new, 2);
	nss2=Allocate_mat2_i(Na_new, 2);
      }
      *N_ali_max=Na_new;
    }

    int dl=lali; if(l1_ini>l2_ini){dl-=l1_ini;}else{dl-=l2_ini;}
    printf("alignment %d vs %d: %d %d -> %d (+%d)\n",
	   n1, n2, l1, l2, lali, dl);
    l1=lali; l2=lali;

    if(0){
      int i;
      printf("a1:"); for(i=0; i<l1; i++)printf(" %d", ali1[i]); printf("\n");
      printf("a2:"); for(i=0; i<l2; i++)printf(" %d", ali2[i]); printf("\n");
    }

    Write_msa(msa, index1, n1, ali1, lali);
    Write_msa(msa, index2, n2, ali2, lali);
  }

  Empty_matrix_i(nc, lmax);
  Empty_matrix_f(d2, lmax);
  Empty_matrix_f(sim_mat, *N_ali_max);
  if(nss1){
    Empty_matrix_i(nss1, *N_ali_max);
    Empty_matrix_i(nss2, *N_ali_max);
  }

  for(int a=0; a<n1; a++)len_msa[index1[a]]=lali;
  for(int a=0; a<n2; a++)len_msa[index2[a]]=lali;

  return(lali);
}

int Align_PC_NW_pair(int *ali_PC, float **d2_in, int **nc_in,
		     struct protein *proti, struct protein *protj)
{
  int l1=proti->len, l2=protj->len;
  int lali=l1+l2, ali1[lali], ali2[lali], vbs=0;
  int **nss1=NULL, **nss2=NULL;
  if(GAP_SS){nss1=Allocate_mat2_i(l1, 2); nss2=Allocate_mat2_i(l2, 2);}
  float pen_pair, penext_pair;

  float **sim_mat=Allocate_mat2_f(l1, l2);
  float **d2=NULL; int **nc=NULL; int comp_nc;
  if(nc_in){ // matrices given from input
    d2=d2_in; nc=nc_in; comp_nc=0;
    pen_pair=penalty; penext_pair=penext;
  }else{
    nc=Allocate_mat2_i(lali, lali); comp_nc=2;
    pen_pair=penalty_ini; penext_pair=penext_ini;
  }

  struct protein *prot_p[2]; prot_p[0]=proti; prot_p[1]=protj;
  Fill_sim_matrix(sim_mat, nss1, nss2, NULL, prot_p, d2, nc, comp_nc,
		  NULL, NULL, NULL, NULL, 1, l1, NULL, 1, l2);

  lali=affinealign_new(sim_mat, nss1, nss2, l1, l2, vbs,
		       pen_pair, penext_pair, ali1, ali2);
  Empty_matrix_f(sim_mat, l1);
  if(nc && nc_in==NULL)Empty_matrix_i(nc, lali);
  if(nss1){
    Empty_matrix_i(nss1, l1);
    Empty_matrix_i(nss2, l2);
  }

  int dl=lali; if(l1>l2){dl-=l1;}else{dl-=l2;}
  if(vbs)printf("pair alignment: %d %d -> %d (+%d)\n", l1, l2, lali, dl);

  for(int i=0; i<l1; i++)ali_PC[i]=-1;
  for(int i=0; i<lali; i++)if(ali1[i]>=0)ali_PC[ali1[i]]=ali2[i];

  return(lali);
}

void Fill_sim_matrix(float **sim_mat, int **nss1, int **nss2,
		     int **msa, struct protein **prot_p,
		     float **d2, int **nc, int comp_nc,
		     int ***Ali_pair, float ****Rot_pair,float ***Shift_pair,
		     int *index1, int n1, int l1,
		     int *index2, int n2, int l2)
{
  if(n1<=0 || n2<=0){
    printf("WARNING in Fill_sim_matrix, n1= %d n2= %d nothing to do\n",
	   n1, n2); return;
  }
  for(int i=0; i<l1; i++)for(int j=0; j<l2; j++)sim_mat[i][j]=0;

  for(int a=0; a<n1; a++){
    int ia, *ires_a;
    if(msa){ia=index1[a]; ires_a=msa[ia];}else{ia=0; ires_a=NULL;}
    struct protein *prot_a=prot_p[ia]; int L_a=prot_a->len;
    for(int b=0; b<n2; b++){
      int ib, *ires_b;
      if(msa){ib=index2[b]; ires_b=msa[ib];}else{ib=1; ires_b=NULL;}
      struct protein *prot_b=prot_p[ib]; int L_b=prot_b->len;
      int i1, i2; if(ia>ib){i1=ib; i2=ia;}else{i1=ia; i2=ib;}
      int *Ali=NULL; float **Rot=NULL, *Shift=NULL;
      if(Ali_pair){Ali=Ali_pair[i2][i1];}
      if(Rot_pair){Rot=Rot_pair[i2][i1]; Shift =Shift_pair[i2][i1];}
      float w[5];
      for(int i=0; i<5; i++){
	w[i]=weight[i]/sqrt(prot_a->max_sim[i]*prot_b->max_sim[i]);
      }
      if(nss1){
	Sum_ss(nss1, prot_a->ss3, ires_a, l1);
	Sum_ss(nss2, prot_b->ss3, ires_b, l2);
      }

      Add_seq_sim(sim_mat, w[1],
		  prot_a->seq_bs, ires_a, l1,
		  prot_b->seq_bs, ires_b, l2);
      if(d2){
        Add_aligned(sim_mat, w[0], L_a, ires_a, l1, L_b, ires_b, l2);
	Add_TM(sim_mat, w[3], d2, Rot, Shift,
	       prot_a, ires_a, l1, prot_b, ires_b, l2);
      }else{
	Add_sec_str(sim_mat, w[2],
		    prot_a->ss3, ires_a, l1,
		    prot_b->ss3, ires_b, l2);
      }
      Add_CO(sim_mat, w[4], nc, Ali, comp_nc,
	     prot_a->Cont_map, L_a, ires_a, l1,
	     prot_b->Cont_map, L_b, ires_b, l2);
    }
  }

  if(NORM_SIM){
    float norm=n1*n2;
    for(int i=0; i<l1; i++)for(int j=0; j<l2; j++)sim_mat[i][j]/=norm;
  }
  return;
}

void Write_msa(int **msa, int *index, int n, int *ali, int lali)
{
  int tmp[lali];
  for(int a=0; a<n; a++){
    int ia=index[a], *msa_a=msa[ia], j; //, last=-1, err=0
    for(j=0; j<lali; j++){
      if(ali[j]>=0){
	tmp[j]=msa_a[ali[j]];
	//if(tmp[j]>=0){if(tmp[j]<=last)err++; last=tmp[j];}
      }else{tmp[j]=-1;}
    }
    /*if(err){
      printf("%d errors in Write_msa seq %d of %d\nali: ", err, ia, n);
      for(j=0; j<lali; j++)printf(" %d", msa_a[j]); printf("\n");
      exit(8);
      }*/
    for(j=0; j<lali; j++){msa_a[j]=tmp[j];}
  }
}

int affinealign_new(float **sim_mat, int **nss1, int **nss2, int  l1, int  l2, 
		    int verbose, float  penalty, float  penext,
		    int *align1, int *align2)
{

  int maxdim = l1; if(l2>maxdim)maxdim=l2;
  
  float **matrix=Allocate_mat2_f(maxdim, maxdim);
  int **dirx=Allocate_mat2_i(maxdim, maxdim);
  int **diry=Allocate_mat2_i(maxdim, maxdim);

  int i, j;
  for(i=0;i<maxdim;i++){
    for(j=0;j<maxdim;j++){
      dirx[i][j] = -1;
      diry[i][j] = -1;
    }
  }

  float gappen[maxdim]; gappen[0]=penalty;
  for(i=1; i<maxdim; i++)gappen[i]=gappen[i-1]+penext;
  
  /* Fill in scores up the right hand side of the matrix   */
  int l1_1=l1-1;
  for(j=0; j<l2; j++){matrix[l1_1][j]=sim_mat[l1_1][j];}

  /* Fill in scores along the bottom row of the matrix  */
  int l2_1=l2-1;        
  for(i=0; i<l1; i++){matrix[i][l2_1]=sim_mat[i][l2_1];}

  i = l1_1;
  j = l2_1;

   /* Move back along the diagonal */
  while(i > 0 && j > 0){
     i--;
     j--;
     
     /* Fill in the scores along this row */
     int j1=j+1;
     for(int i1 = i+1; i1 > 0; i1--){
       float dia = matrix[i1][j1];
       
       /* Find highest score to right of diagonal */
       int rcell=i1+1, gapext = 0; float right, gap_ss=0;
       if(rcell < l1){right=- gappen[1];}else{right=0;}
       for(int k = rcell; k<l1; k++, gapext++){
	 if(nss1){
	   if(nss1[k][0])gap_ss-=nss1[k][0]*SS_sim_score[1][3];
	   if(nss1[k][1])gap_ss-=nss1[k][1]*SS_sim_score[2][3];
	 }
	 float score = matrix[k][j1] - gappen[gapext] - gap_ss;
	 //(penalty + gapext*penext);
	 if(score > right){
	   right = score;
	   rcell = k;
	 }
       }
       
       /* Find highest score below diagonal */
       int dcell=j1+1; float down; gapext=0; gap_ss=0;
       if(dcell < l2){down=- gappen[1];}else{down=0;}
       for(int k = dcell; k<l2; k++, gapext++){
	 if(nss2){
	   if(nss2[k][0])gap_ss-=nss2[k][0]*SS_sim_score[1][3];
	   if(nss2[k][1])gap_ss-=nss2[k][1]*SS_sim_score[2][3];
	 }
	 float score = matrix[i1][k] - gappen[gapext] - gap_ss;
	 if(score > down){
	   down = score;
	   dcell = k;
	 }
       }
         
       /* Set score to best of these */
       float maxoff = right; if(down>maxoff)maxoff=down;
       int ii=i1-1;
       if(dia >= maxoff){
	 matrix[ii][j] = dia;
	 dirx[ii][j] = i1;
	 diry[ii][j] = j1;
       }else if(right > down){
	 matrix[ii][j] = right;
	 dirx[ii][j] = rcell;
	 diry[ii][j] = j1;
       }else{
	 matrix[ii][j] = down;
	 dirx[ii][j] = i1;
	 diry[ii][j] = dcell;
       }
       
       /* Add the score for a match    */  
       matrix[ii][j] +=  sim_mat[ii][j];
       ii=i1;

     } // end i1
     

     /* Fill in the scores in this column */
     int i1=i+1;
     for(int j1 = j+1; j1 >0; j1--){
       float dia  = matrix[i1][j1];
       
       /* Find highest score to right of diagonal */
       int rcell=i1+1, gapext = 0; float right=0, gap_ss=0;
       if(rcell < l1){right=- gappen[1];}
       for(int k = rcell; k<l1; k++, gapext++){
	 if(nss1){
	   if(nss1[k][0])gap_ss-=nss1[k][0]*SS_sim_score[1][3];
	   if(nss1[k][1])gap_ss-=nss1[k][1]*SS_sim_score[2][3];
	 }
	 float score = matrix[k][j1] - gappen[gapext] - gap_ss;
	 if(score > right){
	   right = score;
	   rcell = k;
	 }
       }

       /* Find highest score below diagonal */
       int dcell=j1+1; float down=0; gapext=1; gap_ss=0;
       if(dcell < l2){down= - gappen[1];}
       for(int k = dcell; k<l2; k++, gapext++){
	 if(nss2){
	   if(nss2[k][0])gap_ss-=nss2[k][0]*SS_sim_score[1][3];
	   if(nss2[k][1])gap_ss-=nss2[k][1]*SS_sim_score[2][3];
	 }
	 float score = matrix[i1][k] - gappen[gapext] - gap_ss;
	 if(score > down){
	   down = score;
	   dcell= k;
	 }
       }

       /* Set score to best of these */
       float maxoff = right; if(down>maxoff)maxoff=down;
       int jj=j1-1;
       if(dia >= maxoff){
	 matrix[i][jj] = dia;
	 dirx[i][jj] = i1;
	 diry[i][jj] = j1;
       }else if(right > down){
	 matrix[i][jj] = right;
	 dirx[i][jj] = rcell;
	 diry[i][jj] = j1;
       }else{
	 matrix[i][jj] = down;
	 dirx[i][jj] = i1;
	 diry[i][jj] = dcell;
       }
       
       /* Add the score for a match */
       matrix[i][jj] += sim_mat[i][jj];
       jj=j1;

     } // end j1
  } // end i, j
  
  int align_len;
  float score = TraceBack(matrix, dirx, diry, l1, l2,
			  align1, align2, &align_len);

   if(verbose){
     printf("score= %.4g\n", score);
     printf("Matrix:\n-------\n");
     for(j=l2-1; j>=0;j--){
      printf("j= %d  ",j);
       for(i=l1-1; i>=0; i--){printf(" %.2g",matrix[i][j]);}
       printf("\n");
     }
     printf("Path:\n-----\n");
     for(j=l2-1; j>=0;j--){
       for(i=l1-1; i>=0; i--){
	 printf("(%d,%d) ",dirx[i][j],diry[i][j]);
       }
       printf("\n");
     }
     //exit(8);
   }
    
   Empty_matrix_f(matrix, maxdim);
   Empty_matrix_i(dirx, maxdim);
   Empty_matrix_i(diry, maxdim);

   return(align_len);
}

static float TraceBack(float **matrix, int **dirx, int **diry,
                     int len1, int len2,
                     int *align1, int *align2, int *align_len)
{
   int   i,j,ai,BestI,BestJ,next_x, next_y;

   ai = SearchForBest(matrix, len1, len2, &BestI, &BestJ, align1, align2);
   //printf("starting gap= %d\n", ai);

   /* Now trace back to find the alignment                              */
   i          = BestI;
   j          = BestJ;
   align1[ai] = i;
   align2[ai] = j;
   ai++;

   int len1_1=len1-1, len2_1=len2-1; 
   while(i < len1_1 && j < len2_1){
      next_x = dirx[i][j];
      next_y = diry[i][j];
      i++; j++;
      if((next_x == i) && (next_y == j)){
         /* We are inheriting from the diagonal */
         
      }else if(next_y == j){
         /* We are inheriting from the off-diagonal inserting a gap in
            the y-sequence (seq2) */
         while((i < next_x) && (i < len1_1)){
            align1[ai] = i;
            align2[ai] = -1;
	    i++; ai++;
         }
      }else if(next_x == i){
         /* We are inheriting from the off-diagonal inserting a gap in
            the x-sequence (seq1) */
         while((j < next_y) && (j < len2_1)){
            align1[ai] = -1;
            align2[ai] = j;
	    ai++; j++;
         }
      }else{
         /* Cockup!  */
         fprintf(stderr,"align.c/TraceBack() internal error\n");
	 printf("TraceBack() internal error\n"); exit(8);
      }
      
      align1[ai] = i;
      align2[ai] = j;
      ai++;
   }

   /* If one sequence finished first, fill in the end with insertions   */
   if(i < len1_1){
      for(j=i+1; j<len1; j++){
         align1[ai] = j;
         align2[ai] = -1;
	 ai++;
      }
   }else if(j < len2_1){
      for(i=j+1; i<len2; i++){
	align1[ai] = -1;
	align2[ai] = i;
	ai++;
      }
   }
   
   *align_len = ai;
   
   return(matrix[BestI][BestJ]);
}


static int SearchForBest(float  **matrix, int  len1, int  len2, 
                         int  *BestI, int  *BestJ, int *align1, int *align2)
{
  int   ai=0, i, j, besti=0, bestj=0;
   
   /* Now search the outside of the matrix for the highest scoring cell */
   for(i=1; i<len1; i++){if(matrix[i][0] > matrix[besti][0]) besti = i;}
   for(j=1; j<len2; j++){if(matrix[0][j] > matrix[0][bestj]) bestj = j;}
   if(matrix[besti][0] > matrix[0][bestj]){
      *BestI = besti;
      *BestJ = 0;
      for(i=0; i<*BestI; i++){
         align1[ai] = i;
         align2[ai] = -1; ai++;
      }
   }else{
     *BestI = 0;
     *BestJ = bestj;
     for(j=0; j<*BestJ; j++){
         align1[ai] = -1;
         align2[ai] = j; ai++;
     }
   }
   return(ai);
}

void Set_scores(struct protein *prot)
{
  if(ini_ss){
    float gap=SS_sim_score[0][3];
    for(int i=0; i<3; i++){
      SS_sim_score[i][3]-=gap; 
      SS_sim_score[3][i]-=gap; 
    }
    ini_ss=0;
  }

  prot->max_sim=malloc(5*sizeof(int));
  int L=prot->len, i, error=0;
  float sum_seq=0, sum_ss=0;
  short *seq_bs=prot->seq_bs, *ss3=prot->ss3;
  for(i=0; i<L; i++){
    sum_seq+=blosum62[*seq_bs][*seq_bs];
    sum_ss+=SS_sim_score[*ss3][*ss3];
    if(*seq_bs<0 || *seq_bs >= 22)error++;
    if(*ss3<0 || *ss3 >= 4)error++;
    seq_bs++; ss3++;
  }
  if(error || isnan(sum_seq) || isnan(sum_ss) || sum_seq<=0 || sum_ss<=0){
    printf("ERROR Protein %s L=%d %d errors in seq or secstr, "
	   "<seq_sim>= %.3g <sec_str_sim>= %.3g. Exiting\n",
	   prot->name, prot->len, error,
	   sum_seq/(float)L, sum_ss/(float)L);
    exit(8);
  }

  prot->max_sim[0]=1; // aligned fraction
  prot->max_sim[1]=sum_seq/(float)L; // sequence similarity
  prot->max_sim[2]=sum_ss/(float)L; // secondary structure
  prot->max_sim[3]=1; // TM score
  prot->max_sim[4]=prot->N_cont/(float)L; // contact overlap
}

void Add_aligned(float **sim_mat, float score,
		 int len1, int *ires_1, int l1,
		 int len2, int *ires_2, int l2)
{
  int err=0;
  for(int i1=0; i1<l1; i1++){
    int i; if(ires_1){i=ires_1[i1]; if(i<0)continue;}else{i=i1;}
    if(i>=len1){
      printf("ERROR, site %d >= len= %d (1, %d of %d)\n", i, len1, i1, l1);
      err++;
    }
    for(int i2=0; i2<l2; i2++){
      int j; if(ires_2){j=ires_2[i2]; if(j<0)continue;}else{j=i2;}
      if(j>=len2){
	printf("ERROR, site %d >= len= %d (2, %d of %d)\n", j, len2, i2, l2);
	err++;
      }
      sim_mat[i1][i2]+=score; // Add one if aligned
    }
    if(err){
      printf("ERROR in Add_aligned len1=%d len2=%d, exiting\n", len1, len2);
      int i;
      printf("a1:"); for(i=0; i<l1; i++)printf(" %d", ires_1[i]); printf("\n");
      printf("a2:"); for(i=0; i<l2; i++)printf(" %d", ires_2[i]); printf("\n");
      exit(8);
    }
  }
}

void Add_seq_sim(float **sim_mat, float score,
		 short *seq1, int *ires_1, int l1,
		 short *seq2, int *ires_2, int l2)
{
  for(int i1=0; i1<l1; i1++){
    int i; if(ires_1){i=ires_1[i1]; if(i<0)continue;}else{i=i1;}
    int s1=seq1[i];
    if(s1<0 || s1>=22){
      printf("ERROR in Add_Seq, i=%d seq=%d\n",i, s1); exit(8);
    }
    for(int i2=0; i2<l2; i2++){
      int j; if(ires_2){j=ires_2[i2]; if(j<0)continue;}else{j=i2;}
      int s2=seq2[j];
      if(s2<0 || s2>=22){
	printf("ERROR in Add_Seq, j=%d seq=%d\n", j, s2); exit(8);
      }
      sim_mat[i1][i2]+=score*blosum62[s1][s2]; // blosum score
    }
  }
}

void Add_sec_str(float **sim_mat, float score,
		 short *ss1, int *ires_1, int l1,
		 short *ss2, int *ires_2, int l2)
{
  for(int i1=0; i1<l1; i1++){
    int i, s1; if(ires_1){i=ires_1[i1];}else{i=i1;}
    if(i<0){if(GAP_SS==0)continue; s1=3;}
    else{s1=ss1[i];}
    for(int i2=0; i2<l2; i2++){
      int j, s2; if(ires_2){j=ires_2[i2];}else{j=i2;}
      if(j<0){if(GAP_SS==0)continue; if(i<0)continue; s2=3;}
      else{s2=ss2[j];}
      sim_mat[i1][i2]+=score*SS_sim_score[s1][s2];
    }
  }
}

void Add_TM(float **sim_mat, float score, float **d2,
	    float **Rot, float *Shift,
	    struct protein *prot1, int *ires_1, int l1,
	    struct protein *prot2, int *ires_2, int l2)
{
  float *xca1_store=prot1->xca_rot; int len1=prot1->len;
  float *xca2_store=prot2->xca_rot; int len2=prot1->len;

  int ltar=len1; if(l2>ltar)ltar=len2;
  float d02=TM_coeff*pow(ltar-L_offset, 1./3)-TM_offset; d02*=d02;

  // Rotate structure 2
  if(Shift){ // Rigid body operations are input
    float xca1[3*len1], xca2[3*len2];
    Copy_coord(xca1, xca1_store, 3*len1);
    Copy_coord(xca2, xca2_store, 3*len2);
    Rot_shift(xca2, len2, Rot, Shift);
    All_distances(d2, xca1, len1, xca2, len2);
  }else if(ires_1){ // alignment is input, compute TM score-> rigid body
    if(0){ // Now the coordinates are already rotated
      int ali_ij[len1], i;
      for(i=0; i<len1; i++)ali_ij[i]=-1;
      for(i=0; i<l1; i++){int k=ires_1[i]; if(k>=0)ali_ij[k]=ires_2[i];}
      TM_score(d2, &d02, Rot, Shift, ali_ij, len1, prot1, prot2, 0);
    }
    All_distances(d2, xca1_store, len1, xca2_store, len2);
  } // otherwise, d2 is input for pairwise alignment 

  for(int i1=0; i1<l1; i1++){
    int i; if(ires_1){i=ires_1[i1]; if(i<0)continue;}else{i=i1;}
    for(int i2=0; i2<l2; i2++){
      int j; if(ires_2){j=ires_2[i2]; if(j<0)continue;}else{j=i2;}
      sim_mat[i1][i2]+=score*(1./(1+d2[i][j]/d02)); //exp(-d2[i][j]/d02);
    }
  }

}

void Add_CO(float **sim_mat, float score, int **nc, int *ali, int comp_nc,
	    short **Cont_map1, int len1, int *ires_1, int l1,
	    short **Cont_map2, int len2, int *ires_2, int l2)
{
  int *ali_ij=NULL;
  if(ali){ali_ij=ali;} // pairwise alignment is input
  else if(ires_1){ // multiple alignment is input
    ali_ij=malloc(len1*sizeof(int)); int i;
    for(i=0; i<len1; i++)ali_ij[i]=-1;
    for(i=0; i<l1; i++){int k=ires_1[i]; if(k>=0)ali_ij[k]=ires_2[i];}
  } // otherwise nc is input from pairwise alignment
  if(ali_ij){
    Shared_contacts(nc, ali_ij, Cont_map1, len1, Cont_map2, len2);
    if(ali==NULL)free(ali_ij);
  }else if(comp_nc){
    Shared_contacts_noali(nc, Cont_map1, len1, Cont_map2, len2);
  }
  for(int i1=0; i1<l1; i1++){
    int i; if(ires_1){i=ires_1[i1]; if(i<0)continue;}else{i=i1;}
    for(int i2=0; i2<l2; i2++){
      int j; if(ires_2){j=ires_2[i2]; if(j<0)continue;}else{j=i2;}
      sim_mat[i1][i2]+=score*nc[i][j];
    }
  }
}

void Sum_ss(int **nss, short *ss3, int *ires, int l)
{
  for(int i=0; i<l; i++){
    int k; if(ires){k=ires[i]; if(k<0)continue;}else{k=i;} 
    if(ss3[k]==1){nss[i][0]++;} // helix
    else if(ss3[k]==2){nss[i][1]++;} // strand
  }
}

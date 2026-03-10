#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "protein.h"
#include "tm_score.h"

extern int Align_PC_NW_pair(int *ali_PC, float **d2_in, int **nc_in,
			    struct protein *proti, struct protein *protj);
static float Sequence_identity(int *ali_ij, char *seqi, int li,
			       char *seqj, int lj);
static void Change_ali(int *ali_tmp, int lk, int *ali_ij, int li, int lj);

int Cluster_seqs(int **rep_str, // repr. str. (one per cluster)
		 int **N_all, // total conf per clus
		 int ***rep_conf, int ***rep_index, int **N_conf,
		 // repr. conf (N_conf per clus)
		 int ****all_conf, int ***N_conf_all, // (N_conf per clus)
		 int ****ali_clus, int **L_ali, // Mult. alignments
		 int ****ali_pairs,
		 float ****TM_pairs, // TM scores of cluster i 
		 struct protein *prots, int N_pdb, // input
		 float seq_thr, float str_thr);

int Cluster_seqs(int **rep_str, // repr. str. (one per cluster)
		 int **N_all, // total conf per clus
		 int ***rep_conf, int ***rep_index, int **N_conf,
		 // repr. conf (N_conf per clus)
		 int ****all_conf, int ***N_conf_all, // (N_conf per clus)
		 int ****ali_clus, int **L_ali, // Mult. alignments
		 int ****ali_pairs,
		 float ****TM_pairs, // TM scores of cluster i 
		 struct protein *prots, int N_pdb, // input
		 float seq_thr, float str_thr)
{
  int L_max=0, i, j;
  *ali_pairs=malloc(N_pdb*sizeof(int **));
  for(i=0; i<N_pdb; i++){
    if(i){(*ali_pairs)[i]=malloc(i*sizeof(int *));} //i1>i2
    for(j=0; j<i; j++){(*ali_pairs)[i][j]=NULL;}
  } 
  for(i=0; i<N_pdb; i++){if(prots[i].len>L_max){L_max=prots[i].len;}}
  int N_ali_max=2*L_max;
  int *ali_ij;
    
  int num_all[N_pdb], rep[N_pdb], *conf[N_pdb];
  int change_ali[N_pdb], *ali_tmp[N_pdb];
  for(i=0; i<N_pdb; i++){
    num_all[i]=0; rep[i]=-1; change_ali[i]=0;
    conf[i]=malloc(N_pdb*sizeof(int));
    ali_tmp[i]=malloc(N_ali_max*sizeof(int));
  }

  int N_seq=0;
  printf("Aligning: ");
  for(j=0; j<N_pdb; j++){  // New protein to be classified j>i
    struct protein *prot1=prots+j; int join=0;
    printf(" %d", j); fflush(NULL);
    for(i=0; i<N_seq; i++){ // Already classified proteins i < j
      struct protein *prot2=prots+rep[i];
      if((*ali_pairs)[N_seq][i]==NULL){
	(*ali_pairs)[N_seq][i]=malloc(N_ali_max*sizeof(int));
      }
      ali_ij=(*ali_pairs)[N_seq][i]; // i1=N_seq > i2=i
      Align_PC_NW_pair(ali_ij, NULL, NULL, prot2, prot1);
      float SI=Sequence_identity(ali_ij, prot2->aseq, prot2->len,
				 prot1->aseq, prot1->len);
      if(SI>seq_thr){ // Join to cluster i
	if(0 && prot1->len > prot2->len){ // Change representative
	  rep[i]=j;
	  change_ali[i]=1;
	  for(int a=0; a<num_all[i]; a++){
	    int k=conf[i][a];
	    Change_ali(ali_tmp[k],prots[k].len,ali_ij,prot2->len,prot1->len);
	  }
	  for(int k=0; k<prot1->len; k++){ali_tmp[j][k]=k;} // ident alignment
	}else{
	  for(int k=0; k<prot2->len; k++){ali_tmp[j][k]=ali_ij[k];} 
	}
	conf[i][num_all[i]]=j;
	num_all[i]++; join=1;
	break;
      }
    } // end i
    if(join==0){
      // Start new cluster
      conf[N_seq][0]=j; num_all[N_seq]=1; rep[N_seq]=j;
      for(int k=0; k<prot1->len; k++){ali_tmp[j][k]=k;} // ident alignment
      N_seq++;
    }
  } // end j
  printf("\n%d proteins classified in %d sequence clusters at %.2f identity\n",
	 N_pdb, N_seq, seq_thr);

  // Allocate output
  *N_all=malloc(N_seq*sizeof(int));
  *N_conf=malloc(N_seq*sizeof(int));
  *rep_str=malloc(N_seq*sizeof(int));
  *rep_conf=malloc(N_seq*sizeof(int *));
  *rep_index=malloc(N_seq*sizeof(int *));
  *all_conf=malloc(N_seq*sizeof(int **));
  *N_conf_all=malloc(N_seq*sizeof(int *));

  *ali_clus=malloc(N_seq*sizeof(int **));
  *L_ali=malloc(N_seq*sizeof(int));
  *TM_pairs=malloc(N_seq*sizeof(float **));

  for(i=0; i<N_seq; i++){
    int n=num_all[i];
    int L_max=-1;
 
    printf("Cluster %d, n=%d\n", i, n);
    (*N_all)[i]=n;
    (*ali_clus)[i]=malloc(n*sizeof(int *));
    int kr=rep[i];
    struct protein *protr=prots+kr, *prot_p[n];
    (*L_ali)[i]=protr->len;
    for(j=0; j<n; j++){
      int k=conf[i][j];
      if(k<0 || k>=N_pdb){
	printf("ERROR in Cluster_seq, wrong conformation index %d (%d)\n",
	       k, N_pdb); exit(8);
      }
      prot_p[j]=prots+k;
      int L=prots[k].len; if(L>L_max){L_max=L;}
    }

    //printf("Representative structure: %d change: %d\n", kr, change_ali[i]);
    // Do multiple alignment with current rep_str
    for(j=0; j<n; j++){
      int k=conf[i][j], *ali_k=ali_tmp[k], l;
      (*ali_clus)[i][j]=malloc(L_max*sizeof(int));
      if(change_ali[i] && k!=kr){
	Align_PC_NW_pair(ali_k, NULL, NULL, protr, prots+k);
      }
      if(k!=kr){
	for(l=0; l<protr->len; l++){(*ali_clus)[i][j][l]=ali_k[l];}
      }else{
	for(l=0; l<protr->len; l++){(*ali_clus)[i][j][l]=l;}
      }
    }
    // Determine repr str as most central str
    (*TM_pairs)[i]=malloc(n*sizeof(float *));
    float **TM_i= (*TM_pairs)[i], sum_TM[n];
    for(j=0; j<n; j++){
      TM_i[j]=malloc(n*sizeof(float));
      TM_i[j][j]=1.;
    }
    int r=0;
    if(n>1){
      float TMs=TM_score_mult(TM_i,(*ali_clus)[i],(*L_ali)[i],n,prot_p,0);
      for(j=0; j<n; j++){
	sum_TM[j]=0;
	for(int k=0; k<j; k++){
	  sum_TM[j]+=TM_i[j][k];
	  sum_TM[k]+=TM_i[j][k];
	}
      }
      for(j=1; j<n; j++){
	if(sum_TM[j]>sum_TM[r]){r=j;}
      }
      kr=conf[i][r];
      protr=prots+kr;
      (*L_ali)[i]=protr->len;
    }
    (*rep_str)[i]=kr;
    //printf("Representative structure: %d (%d) previous: %d\n", kr, r, rep[i]);

    // Determine multiple alignment vs representative str
    if(kr!=rep[i]){
      for(j=0; j<n; j++){
	int k=conf[i][j], l;
	if(k!=kr){
	  Align_PC_NW_pair((*ali_clus)[i][j],NULL,NULL,protr,prots+k);
	}else{
	  for(l=0; l<protr->len; l++){(*ali_clus)[i][j][l]=l;}
	}
	for(l=protr->len; l<L_max; l++){(*ali_clus)[i][j][l]=-1;}
      }
    }

    // Determine conformations
    int nc=0, nconf[n];
    if(n==1){nc=1; nconf[0]=0;}
    else{
      for(j=0; j<n; j++){
	int k=0, join=0;
	for(k=0; k<j; k++){
	  if(TM_i[j][k]>str_thr){ // Join j to cluster k
	    nconf[j]=nconf[k]; join=1; break;
	  }
	}
	if(join==0){nconf[j]=nc; nc++;}
      }
      // Find representatives
      for(j=0; j<n; j++){
	sum_TM[j]=0;
	for(int j2=0; j2<j; j2++){
	  if(nconf[j2]==nconf[j]){
	    sum_TM[j] +=TM_i[j][j2];
	    sum_TM[j2]+=TM_i[j][j2];
	  }
	}
      }
    }

    printf("%d conformations\n", nc);
    printf("Representative structure: %d (%d) previous: %d\n", kr, r, rep[i]);
    (*N_conf)[i]=nc;
    (*all_conf)[i]=malloc(nc*sizeof(int *));
    (*N_conf_all)[i]=malloc(nc*sizeof(int));
    (*rep_conf)[i]=malloc(nc*sizeof(int));
    (*rep_index)[i]=malloc(nc*sizeof(int));
    for(int k=0; k<nc; k++){
      (*all_conf)[i][k]=malloc(n*sizeof(int));
      int r=-1, m=0;
      for(j=0; j<n; j++){
	if(nconf[j]!=k){continue;}
	(*all_conf)[i][k][m]=conf[i][j]; m++;
	if(r<0 || sum_TM[j]>sum_TM[r]){r=j;}
      }
      if(r<0){
	printf("ERROR, conf %d of %d has no structure\n", k, nc); exit(8);
      }
      //printf("conf %d has %d structure, rep=%d (%d)\n", k, m, r, conf[i][r]); 
      (*rep_conf)[i][k]=conf[i][r];
      (*rep_index)[i][k]=r;
      (*N_conf_all)[i][k]=m;
    }
    
  }

  int change=0;
  for(int i1=0; i1<N_seq; j++){
    struct protein *prot1=prots+(*rep_str)[i1];
    for(int i2=0; i2<i1; i2++){
      if((*rep_str)[i1]!=rep[i1] || (*rep_str)[i2]!=rep[i2]){
	change++; printf("change= %d\n", change);
	// Change pairwise alignments
	ali_ij=(*ali_pairs)[i1][i2]; // i1 >i2
	Align_PC_NW_pair(ali_ij,NULL,NULL,prots+(*rep_str)[i2],prot1);
      }
    }
  }

  // Clean
  for(i=0; i<N_pdb; i++){
    free(conf[i]);
    free(ali_tmp[i]);
  }
  return(N_seq);
  }

void Change_ali(int *ali_tmp, int lk, int *ali_ij, int li, int lj)
{
  // change ali_tmp aligned to protein li to protein lj aligned with ali_ij
  int tmp[lj]; for(int i=0; i<lj; i++){tmp[i]=-1;}

  for(int i=0; i<li; i++){
    int a=ali_ij[i]; if(a>=0)tmp[a]=ali_tmp[i];
  }
  for(int i=0; i<lj; i++){ali_tmp[i]=tmp[i];}
  return;
}

float Sequence_identity(int *ali_ij, char *seqi, int li, char *seqj, int lj)
{
  int id=0, lmin=li; if(lj<lmin){lmin=lj;}
  for(int k=0; k< li; k++){
    int k2=ali_ij[k];
    if(k2>=0 && seqi[k]==seqj[k2]){id++;}
  }
  return((float)id/lmin);
}

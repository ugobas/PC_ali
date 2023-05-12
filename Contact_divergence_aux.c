#include "Contact_divergence_aux.h"
#include "normalization.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "NeedlemanWunsch.h"
#include "allocate.h"
int VERBOSE=1;  // Verbose: 


static int Equivalent(char aa1, char aa2);
static char *Remove_gaps(char *gapped_seq, int N_ali, int *N);
void Max_shared(int *neighb1, float *score1, int n1,
		int *neighb2, float *score2, int n2,
		int **nc, int *ali);


int Align_seq(int *Prot_ali, int N_ali, char *Prot_seq, char *PDB_seq, int n2)
{
  // Prot_ali: position l in the PDB aligned to i in seq
  int l=0, i, Lali=0, lPDB=n2, gap=0, noali=0;
  for(i=0; i<N_ali; i++){
    Prot_ali[i]=-1;
    if(Prot_seq[i]!='-'){
      if(Equivalent(Prot_seq[i],PDB_seq[l])){
	Prot_ali[i]=l; l++; Lali++;
      }else if(l==0){
	while((l<lPDB)&&(Equivalent(Prot_seq[i],PDB_seq[l]==0))){gap++; l++;}
      }else{
	noali++;
      }
    }
  }
  if(VERBOSE)printf("%d residues perfectly aligned with %d gaps in PDB\n",
		    Lali, gap);
  if(noali){
    printf("%d gaps in the middle of PDB sequence, leaving\n",noali);
    return(-1); //goto align;
  }
  return(Lali);

  // Perform sequence alignment
  Lali=0;
  // Sequence alignment parameters
  int IDE=1;  // Use identity to score alignment
  int GAP=7;  // Gap opening penalty
  int n1, nali;
  char *seq_nogap=Remove_gaps(Prot_seq, N_ali, &n1);
  printf("Aligning %d res in ali and %d in PDB\n",n1,n2);
  char ali1[n1+n2], ali2[n1+n2];
  int verbose=0;
  if(alignNW(seq_nogap,n1,PDB_seq,n2,verbose,IDE,GAP,ali1,ali2,&nali)==0){
    printf("ERROR, alignment failed\n"); return(-1);
  }
  int mismatch=0, nomatch[n1], lnomatch[n1];
  int l1=-1, l2=0; gap=0;
  int ali[n1];
  for(i=0; i<nali; i++){
    if(ali1[i]!='-'){
      l1++;
      if(ali2[i]!='-'){
	ali[l1]=l2; Lali++;
	if(Equivalent(ali1[i],ali2[i])==0){
	  printf("Mismatch: %c-%c i=%d\n", ali1[i], ali2[i], i);
	  nomatch[mismatch]=i; lnomatch[mismatch]=l1;
	  mismatch++;
	}
      }else{
	ali[l1]=-1;  gap++;
      }
    }
    if(ali2[i]!='-')l2++;
  }
  l1=0;
  for(i=0; i<N_ali; i++){
    if(Prot_seq[i]=='-'){
      Prot_ali[i]=-1;
    }else{
      Prot_ali[i]=ali[l1]; l1++;
    }
  }
  if(VERBOSE){
    for(i=0; i<N_ali; i++){printf("%c", Prot_seq[i]);} printf("\n");
    for(i=0; i<N_ali; i++){
      int k=Prot_ali[i];
      if(k>=0){printf("%c", PDB_seq[k]);}else{printf("-");}
    }
    printf("\n");
    printf("%d residues aligned with %d gaps and %d mismatches: ",
	   Lali, gap, mismatch);
    for(i=0; i<mismatch; i++){
      int k=nomatch[i];
      printf(" %c%d%c", ali1[k], lnomatch[i], ali2[k]);
    }
    printf("\n");
  }
  free(ali1); free(ali2); free(seq_nogap);
  return(Lali);
}

static char *Remove_gaps(char *gapped_seq, int N_ali, int *N){
  *N=0;
  char *seq=malloc(N_ali*sizeof(char));
  for(int i=0; i<N_ali; i++){
    if(gapped_seq[i]!='-'){
      seq[*N]=gapped_seq[i]; (*N)++;
    }
  }
  return(seq);
}

static int Equivalent(char aa1, char aa2){
  if((aa1==aa2)||(aa1=='X')||(aa2=='X'))return(1);
  return(0);
}

float Seq_identity(char *seq1, char *seq2, int N_ali){
  int i, id=0, l1=0, l2=0;
  for(i=0; i<N_ali; i++){
    if(seq1[i]!='-'){
      l1++; if(seq1[i]==seq2[i])id++;
    }
    if(seq2[i]!='-')l2++; 
  }
  float lnorm;
  if(NORM==0){if(l1<l2){lnorm=l1;}else{lnorm=l2;}}
  else if(NORM==1){if(l1>l2){lnorm=l1;}else{lnorm=l2;}}
  else{lnorm=sqrt(l1*l2);}
  return(id/lnorm);
}

void Align_CO(int *ali_CO, int *ali, int **nc,
	      short **Cont_map1, int n1,
	      short **Cont_map2, int n2)
{
  int neighb1[n1], neighb2[n2], ali_tmp[n1];
  float thr=4; // Minimum number of contacts for alignment
  float score1[n1], score2[n2];

  Shared_contacts(nc, ali, Cont_map1, n1, Cont_map2, n2);
  Max_shared(neighb1, score1, n1, neighb2, score2, n2, nc, ali);
  Align_neighbors(ali_tmp, ali, thr, neighb1, score1, n1, neighb2,score2,n2);
  Shared_contacts(nc, ali_tmp, Cont_map1, n1, Cont_map2, n2);
  Max_shared(neighb1, score1, n1, neighb2, score2, n2, nc, ali);
  Align_neighbors(ali_CO, ali_tmp, thr, neighb1,score1,n1,neighb2,score2,n2);
}

void Shared_contacts(int **nc, int *ali,
			 short **Cont_mat_1, int n1,
			 short **Cont_mat_2, int n2)
{
  int i1, i2;
  for(i1=0; i1<n1; i1++)for(i2=0; i2<n2; i2++)nc[i1][i2]=0;

  for(i1=0; i1<n1; i1++){
    for(i2=0; i2<n2; i2++){
      short *j1=&Cont_mat_1[i1][0];   
      while(*j1>=0){
	short *j2=&Cont_mat_2[i2][0];
	while(*j2>=0){
	  if(ali[i1]==i2)nc[*j1][*j2]++;
	  if(ali[i1]==*j2)nc[*j1][i2]++;
	  if(ali[*j1]==i2)nc[i1][*j2]++;
	  if(ali[*j1]==*j2)nc[i1][i2]++;
	  j2++;
	}
	j1++;
      }
    }
  }
}

void Shared_contacts_2(int **nc, int *ali,
		       short **Cont_list_1, int n1,
		       short **Cont_list_2, int n2)
{
  int i1, i2;
  for(i1=0; i1<n1; i1++)for(i2=0; i2<n2; i2++)nc[i1][i2]=0;
  
  int **Cont_mat_2=Allocate_mat2_i(n2, n2);
  for(i2=0; i2<n2; i2++){
     short *j2=&Cont_list_2[i2][0];
     while(*j2>=0){Cont_mat_2[i2][*j2]=1; j2++;}
  }

  for(i1=0; i1<n1; i1++){
    i2=ali[i1]; if(i2<0 || i2>=n2)continue;
    short *j1=&Cont_list_1[i1][0];
    while(*j1>=0){
      int j2=ali[*j1];
      if(j2<0 || j2>=n2){j1++; continue;}
      if(Cont_mat_2[i2][j2]){
	nc[i1][i2]++;
	nc[*j1][j2]++;
      }
      j1++;
    }
  }
  Empty_matrix_i(Cont_mat_2, n2);

}

void Shared_contacts_noali(int **nc,
			   short **Cont_list_1, int n1,
			   short **Cont_list_2, int n2)
{
  int i1, i2;
  for(i1=0; i1<n1; i1++)for(i2=0; i2<n2; i2++)nc[i1][i2]=0;

  int **Cont_mat_2=Allocate_mat2_i(n2, n2);
  for(i2=0; i2<n2; i2++){
     short *j2=&Cont_list_2[i2][0];
     while(*j2>=0){Cont_mat_2[i2][*j2]=1; j2++;}
  }

  for(i1=0; i1<n1; i1++){
    for(i2=0; i2<n2; i2++){
      int i12=i2-i1;
      short *j1=&Cont_list_1[i1][0];
      while(*j1>=0){
	int j2=i12+*j1;
	if(Cont_list_2[i2][j2]){
	  nc[i1][i2]++;
	  nc[*j1][j2]++;
	}
	j1++;
      }
    }
  }
  Empty_matrix_i(Cont_mat_2, n2);

}

void Max_shared(int *neighb1, float *score1, int n1,
		int *neighb2, float *score2, int n2,
		int **nc, int *ali)
{
  // nthr: Minimum number of contacts to establish correspondence
  int i1, i2;
  for(i1=0; i1<n1; i1++){
    int imax=-1, nmax=0, a=ali[i1];
    for(i2=0; i2<n2; i2++){
      if(nc[i1][i2]>nmax){
	nmax=nc[i1][i2]; imax=i2;
      }else if(nc[i1][i2]==nmax && imax>=0 && a>=0 && abs(i2-a)<abs(imax-a)){
	nmax=nc[i1][i2]; imax=i2;
      }
    }
    score1[i1]=nmax;
    neighb1[i1]=imax;
  }

  int ali2[n2]; Invert_ali(ali2, n2, ali, n1);
  for(i2=0; i2<n2; i2++){
    int imax=-1, nmax=0, a=ali2[i2];
    for(i1=0; i1<n1; i1++){
      if(nc[i1][i2]>nmax){
	nmax=nc[i1][i2]; imax=i1;
      }else if(nc[i1][i2]==nmax && imax>=0 && a>=0 && abs(i1-a)<abs(imax-a)){
	nmax=nc[i1][i2]; imax=i1;
      }
    }
    score2[i2]=nmax;
    neighb2[i2]=imax;
  }

}

void Invert_ali(int *ali2, int n2, int *ali, int n1){
  for(int i=0; i<n2; i++)ali2[i]=-1;
  for(int i=0; i<n1; i++)if(ali[i]>=0)ali2[ali[i]]=i;
}

float Contact_overlap(int *ali,
		      short **Cont_mat_1, int N1,
		      short **Cont_mat_2, int N2)
{
  int nc1=0, nc2=0, nc12=0, ncali=0, i;
  for(i=0; i<N1; i++){
    short *j1=&Cont_mat_1[i][0];
    while(*j1>=0){nc1++; j1++;}
  }
  for(i=0; i<N2; i++){
    short *j1=&Cont_mat_2[i][0];
    while(*j1>=0){nc2++; j1++;}
  }

  int ali2[N2]; Invert_ali(ali2, N2, ali, N1);
  for(int i1=0; i1<N1; i1++){
    int i2=ali[i1]; if(i2<0)continue;
    short *j1=&Cont_mat_1[i1][0];   
    short *j2=&Cont_mat_2[i2][0];
    while(*j1>=0){
      int j2ali=ali[*j1];
      if(j2ali < 0){j1++; continue;}
      while(*j2>=0){
	if(*j2==j2ali){
	  nc12++;
	  j2++; break;
	}
	j2++;
      } // end j2
      j1++;
    } // end j1
    while(*j2>=0){
      int j1ali=ali2[*j2];
      if(j1ali<0){j2++; continue;}
      ncali++;
      j2++;
    } // end j2
  } // end i1

  float nc; // Normalization
  if(NORMA){nc=ncali;}
  else if(NORM==0){if(nc1<nc2){nc=nc1;}else{nc=nc2;}}
  else if(NORM==1){if(nc1>nc2){nc=nc1;}else{nc=nc2;}}
  else{nc=sqrt(nc1*nc2);}
  float q=nc12; if(nc){q/=nc;};
  return(q);
}

float Test_contacts(int *ali, int *id_3D, float c_ave,
		    int *ali_cont_ct,  int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid,int *id_cont_aaid, 
		    short **Cont_mat_1, char *seq1, int *ncont1, int N1,
		    short **Cont_mat_2, char *seq2, int *ncont2, int N2)
{
  short i, *j1, *j2, i2; int ncali=0;

  for(i=0; i<2; i++){ali_cont_ct[i]=0; id_cont_ct[i]=0;}
  for(i=0; i<=2; i++){
    ali_cont_sup[i]=0; id_cont_sup[i]=0;
    ali_cont_aaid[i]=0; id_cont_aaid[i]=0;
  }

  float c4_ave=pow(c_ave, 4);

  int ali2[N2]; Invert_ali(ali2, N2, ali, N1);
  int nc1=0, nc2=0, nc12=0, kc;
  for(i=0; i<N1; i++){
    j1=&Cont_mat_1[i][0];
    while(*j1>=0){nc1++; j1++;}
  }
  for(i=0; i<N2; i++){
    j1=&Cont_mat_2[i][0];
    while(*j1>=0){nc2++; j1++;}
  }

  for(int i1=0; i1<N1; i1++){
    i2=ali[i1]; if(i2<0)continue;
    int ni3D=0, n3D=0, niid=0, nid=0; float ci=0;
    ci=ncont1[i1]*ncont2[i2];
    if(id_3D[i1]>0){ni3D=1;}else{ni3D=0;}
    if(seq1){if(seq1[i1]==seq2[i2]){niid=1;}else{niid=0;}}
 
    j1=&Cont_mat_1[i1][0];   
    j2=&Cont_mat_2[i2][0];
    while(*j1>=0){
      int j2ali=ali[*j1];
      if(j2ali < 0){j1++; continue;}
      if(ali_cont_ct){
	float c4=ci*ncont1[*j1]*ncont2[j2ali];
	if(c4<=c4_ave){kc=0;}else{kc=1;}
	ali_cont_ct[kc]++;
      }
      if(id_3D[*j1]>0){n3D=ni3D+1;}else{n3D=ni3D;}
      ali_cont_sup[n3D]++;
      if(seq1){
	if(seq1[*j1]==seq2[j2ali]){nid=niid+1;}else{nid=niid;}
	ali_cont_aaid[nid]++;
      }
      while(*j2>=0){
	if(*j2==j2ali){
	  nc12++;
	  id_cont_ct[kc]++;
	  id_cont_sup[n3D]++;
	  id_cont_aaid[nid]++;
	  j2++; break;
	}else if(ali2[*j2]>=0){
	  int j1ali=ali2[*j2];
	  float c4=ci*ncont1[j1ali]*ncont2[*j2];
	  if(c4<=c4_ave){ali_cont_ct[0]++;}
	  else{ali_cont_ct[1]++;}
	  if(id_3D[j1ali]>0){ali_cont_sup[ni3D+1]++;}
	  else{ali_cont_sup[ni3D]++;}
	  if(seq1){
	    if(seq1[j1ali]==seq2[*j2]){ali_cont_aaid[niid+1]++;}
	    else{ali_cont_aaid[niid]++;}
	  }
	}
	j2++;
      } // end j2
      j1++;
    } // end j1
    while(*j2>=0){
      int j1ali=ali2[*j2];
      if(j1ali<0){j2++; continue;}
      ncali++;
      float c4=ci*ncont1[j1ali]*ncont2[*j2];
      if(c4<=c4_ave){ali_cont_ct[0]++;}
      else{ali_cont_ct[1]++;}
      if(id_3D[j1ali]>0){ali_cont_sup[ni3D+1]++;}
      else{ali_cont_sup[ni3D]++;}
      if(seq1){
	if(seq1[j1ali]==seq2[*j2]){ali_cont_aaid[niid+1]++;}
	else{ali_cont_aaid[niid]++;}
      }
      j2++;
    } // end j2
  } // end i1
  
  float nc; // Normalization
  if(NORMA){nc=ncali;}
  else if(NORM==0){if(nc1<nc2){nc=nc1;}else{nc=nc2;}}
  else if(NORM==1){if(nc1>nc2){nc=nc1;}else{nc=nc2;}}
  else{nc=sqrt(nc1*nc2);}

  float q=nc12; if(nc)q/=nc;
  return(q);
}


void Align_neighbors(int *ali_new, int *ali, float thr,
		     int *neighb1, float *score1, int n1,
		     int *neighb2, float *score2, int n2)
{
  int neigh_2[n1]; Invert_ali(neigh_2, n1, neighb2, n2);
  int l1=0, l2=-1; // left sites
  /*float sco_2[n1]; for(i=0; i<n1; i++)sco_2[i]=-1;
    for(i=0; i<n1; i++)if(neigh_2[i]>=0)sco_2[i]=score2[neigh_2[i]]; */

  for(int i1=0; i1<n1; i1++){
    int i2=ali[i1]; 
    if(i2 >=0 && i2>l2 && i2<n2 && (i2==neighb1[i1] || neigh_2[i1]==i2)){
      ali_new[i1]=i2;
      
      for(int j1=l1; j1<i1; j1++){
	int j2=neighb1[j1]; if(j2<=l2 || j2>=i2)j2=-1;
	int j3=neigh_2[j1]; if(j3<=l2 || j3>=i2)j3=-1;
	int ja=ali[j1];     if(ja<=l2 || ja>=i2)ja=-1;
	/*float sc2=score1[j1], sc3=sco_2[j1];
	if(j2<0 || (j3>=0 && sc3>=thr &&(sc2<thr||j3<j2))){j2=j3; sc2=sc3;}
	if(ja<0 || (j2>=0 && (sc2>=thr || j2<ja))){ja=j2;}
	*/
	if(j2<0 || (j3>=0 && j3<j2)){j2=j3;}
	if(ja<0 || (j2>=0 && j2<ja)){ja=j2;}
	if(ja>=n2)ja=-1;
	ali_new[j1]=ja;
	if(ja>=0)l2=ja;
      }
      l1=i1+1; l2=i2;
    }else{
      ali_new[i1]=-1;
    }
  }
  for(int j1=l1; j1<n1; j1++){
    int j2=neighb1[j1]; if(j2<=l2 || j2>=n2)j2=-1;
    int j3=neigh_2[j1]; if(j3<=l2 || j3>=n2)j3=-1; 
    int ja=ali[j1];     if(ja<=l2 || ja>=n2)ja=-1;
    /*sc2=score1[j1]; sc3=sco_2[j1];
    if(j2<0 || (j3>=0 && sc3>=thr &&(sc2<thr||j3<j2))){j2=j3; sc2=sc3;}
    if(ja<0 || (j2>=0 && (sc2>=thr || j2<ja))){ja=j2;}
    */
    if(j2<0 || (j3>=0 && j3<j2)){j2=j3;}
    if(ja<0 || (j2>=0 && j2<ja)){ja=j2;}
    if(ja>=n2)ja=-1;
    ali_new[j1]=ja;
    if(ja>=0)l2=ja;
  }

  // int *a=ali_new; for(i=0; i<n1; i++){if((*a)>=n2)*a=-1; a++;} 
}

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

int Pair_ali(int *ali, int N_ali, int *ali1, int *ali2)
{
  // Store Alignment of sequences 1 and 2
  int nali=0;
  for(int i=0; i<N_ali; i++)ali[i]=-1;
  for(int i=0; i<N_ali; i++){
    if(ali1[i]>=0){ali[ali1[i]]=ali2[i]; if(ali2[i]>=0)nali++;}
  }
  return(nali);
}

float Seqid(int *ali_12, int *id_aa,
	    char *seq1, int N1,
	    char *seq2, int N2)
{
  int id=0, nali=0;
  if(id_aa)for(int i1=0; i1<N1; i1++)id_aa[i1]=0;
  for(int i1=0; i1<N1; i1++){
    int i2=ali_12[i1];
    if(i2>=0){
      if(i2>=N2){
	printf("ERROR in Seqid, i2=%d but l2= %d\n", i2, N2); return(-1);
      }
      nali++;
      if(seq1[i1]==seq2[i2]){id++; if(id_aa)id_aa[i1]=1;}
    }
  }

  float lnorm; // Normalization
  if(NORMA){lnorm=nali;}
  else if(NORM==0){if(N1<N2){lnorm=N1;}else{lnorm=N2;}} // Min
  else if(NORM==1){if(N1>N2){lnorm=N1;}else{lnorm=N2;}} // Max
  else{lnorm=sqrt(N1*N2);} //
  return((float)id/lnorm);
}

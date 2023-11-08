/* TM score = 1/L sum_i 1/[1+(di/d0)^2]
   d0=TM_coeff*(L-L_offset)^(1/3)-TM_offset
*/
#define TM_coeff 1.24
#define TM_offset 1.8
#define L_offset 15
//#include "normalization.h"
#include "McLachlan.h"
#include "protein.h"
#include "tm_score.h"
#include "Contact_divergence_aux.h"
#include "allocate.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Parameters
int N_SUP_AVE=2; // Iterartive construction of ave conformation 
int nali_min=16; // Minimal number of aligned residues
int nfail_max=3; // Max number of decrease of TM score

float Compute_TM_mult(double *tm_k, float **TM_all, float **xca,
		      int **msa, int nprot, int n, float d02, int *len);
int Superimpose_mult(float *w, float **xca, int *len, int *mprot,
		     int i_ref, int **msa, int nprot, int n);
void Add_coord(float *xsum, int *sum, int nn, float *xca_i, int *msa_i,
	       int *iali);
void Remove_coord(float *xsum, int *sum, int nn, float *xca_i, int *msa_i,
		  int *iali);

int Set_weights(int *w1, int *w2, float *w_ali, int *ali, int n1, int n2);
int Set_coord_ali(float *xca1, float *xca1_store, 
		  float *xca2, float *xca2_store, int *ali, int n1);

int Set_coord_cm(float *xca1, double *cm, float *xca1_store, int *w, int n1);

float dist2(float *x1, float *x2);
void Rotate_vector(float *x, float **rot);
void Closest_neighbor(int *neighb1, float *dmin1, float *score, float d0,
		      float **d2, int n1, int n2, int order);
void All_distances(float **d2, float *x1, int n1, float *x2, int n2);

float TM_score(float **d2, float *d02, float **Rot_out, float *Shift_out,
	       int *ali, float ltar,
	       struct protein *prot1, struct protein *prot2, int verbose)
{
  if(ltar <= L_offset)return(-1);
  *d02=TM_coeff*pow(ltar-L_offset, 1./3)-TM_offset;
  if(*d02 <= 0)return(-1);
  if(verbose)printf("TM score, d0= %.3f n=%3.0f\n", *d02, ltar);
  (*d02)*=(*d02);

  // Copy coordinates of aligned residues
  float *xca1_store=prot1->xca_rot; int n1=prot1->len;
  float *xca2_store=prot2->xca_rot; int n2=prot2->len;
  float xca1[3*n1], xca2[3*n2];
  int n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
  if(n<4)return(0);

  double TM_max=0; int i, ncart=3*n;
  float w_TM[n], w[n], dd[n];
  for(i=0; i<n; i++)w_TM[i]=1;

  // Optimizing rotations
  for(int L_ini=n; L_ini>=(n/4); L_ini*=0.5){
    if(L_ini<nali_min)break;
    for(int ini=0; ini<n; ini+=L_ini){
      int end=ini+L_ini; if(end>n)break;
      // initial superimposed fragment
      for(i=0; i<n; i++){
	if((i>=ini)&&(i<end)){w[i]=1;}else{w[i]=0;}
      }
      int nan=0;
      for(i=0; i<ncart; i++)if(isnan(xca1[i])||isnan(xca2[i])){nan=1; break;}
      if(nan)Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);

      float TM_max_tmp=0; int nfail=0, align;
      for(int iter=0; iter<n; iter++){
	align=0; for(i=0; i<n; i++)if(w[i])align++;
	if(align<nali_min)break;
	rmsd_mclachlan_f(NULL, NULL, xca1, xca2, w, n);
	double TM=0; float dd_ave=0; int norm=0;
	for(i=0; i<n; i++){
	  int k=3*i; dd[i]=dist2(xca1+k, xca2+k);
	  if(isnan(dd[i]))continue;
	  TM+=1./(1+dd[i]/(*d02));
	  if(dd[i]<(*d02)){dd_ave+=dd[i]; norm++;}
	}
	if(TM>TM_max_tmp){
	  TM_max_tmp=TM; nfail=0;
	  if(TM>TM_max){TM_max=TM; for(i=0; i<n; i++)w_TM[i]=w[i];}
	}else{
	  nfail++; if(nfail>nfail_max)break;
	}
	// Superimpose only if di<thr
	if(norm){dd_ave/=norm;}
	float dd_thr=((*d02)+iter*dd_ave)/(iter+1);
	for(i=0; i<n; i++)if(dd[i]<dd_thr){w[i]=1;}else{w[i]=0;}
      } // end iter

      if(verbose){
	printf("remove d<d0 L=%3d ini=%3d align=%d ltar=%.0f TM= %.2g\n",
	       L_ini, ini, align, ltar, TM_max_tmp/ltar);
      }
    }
  }
  if(verbose){ //
    printf("TM= %.1f d0=%.1f ltar=%.0f\n", TM_max, sqrt(*d02), ltar);
  }

  // Superimpose all residues, including not aligned
  if(Rot_out || d2){
    n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
    float shift[3], *rot[3]; int j;
    for(i=0; i<3; i++)rot[i]=malloc(3*sizeof(float));
    rmsd_mclachlan_f(rot, shift, xca1, xca2, w_TM, n);
    if(Shift_out)for(i=0; i<3; i++)Shift_out[i]=shift[i];
    if(Rot_out)for(i=0; i<3; i++)for(int j=0; j<3; j++)Rot_out[i][j]=rot[i][j];
    //Rot_shift(xca2_store, n2, rot, shift);
    if(d2){
      Copy_coord(xca2, xca2_store, 3*n2);
      Rot_shift(xca2, n2, rot, shift);
      All_distances(d2, xca1_store, n1, xca2, n2);
    }
    for(i=0; i<3; i++)free(rot[i]);
    //Empty_matrix_f(rot, 3);
  }
  if(TM_max>ltar){
    printf("ERROR, TM score = %.1f/%.1f > 1\n", TM_max,ltar);
    exit(8);
  }

  return(TM_max/ltar);
}

float TM_score_mult(float **TM_all, int **msa, int N_ali, int nprot,
		    struct protein **prot_p, int verbose)
{

  float *xca[nprot];
  int len[nprot], ncart[nprot], i, j, k;
  for(i=0; i<nprot; i++){
    len[i]=prot_p[i]->nca;
    ncart[i]=3*len[i];
    xca[i]=malloc(ncart[i]*sizeof(float));
    Copy_coord(xca[i], prot_p[i]->xca_rot, ncart[i]);
  }

  // Choose first aligned structure as longest structure
  int i_ref=0, lmax=len[0], lmin=lmax;
  for(i=1; i<nprot; i++){
    if(len[i]>lmax){lmax=len[i]; i_ref=i;}
    else if(len[i]<lmin){lmin=len[i];}
  }

  // Set parameters of TM score
  if(lmin <= L_offset)return(-1);
  float d02=TM_coeff*pow(lmin-L_offset, 1./3)-TM_offset; //lmin
  if(d02 <= 0)return(-1);
  d02*=d02;

  // mprot[k]= number of proteins aligned at column k
  int n=N_ali, n2=0, mprot[n], m2=nprot/2;
  float tm_thr[n], w_opt[n], w[n];
  for(k=0; k<n; k++){
    int m=0; for(i=0; i<nprot; i++)if(msa[i][k]>=0)m++;
    if(m>1)n2++;
    mprot[k]=m; w_opt[k]=1;
    tm_thr[k]=m*(m-1)/2;
    if(m>=m2){tm_thr[k]*=0.5;}
    else{tm_thr[k]*=0.75;}
  }

  printf("Performing multiple superimposition of %d aligned structures, "
	 "reference structure: %d d0= %.1f lmin=%d nali=%d col. w pairs: %d\n",
	 nprot, i_ref, sqrt(d02), lmin, n, n2);

  // Optimizing superimposition
  float *rot[3]; for(i=0; i<3; i++)rot[i]=malloc(3*sizeof(float));
  double TM_max=0, tm_k[n];
  for(int L_ini=n; L_ini>=(n/4); L_ini/=2){
    if(L_ini<nali_min)break;
    for(int ini=0; ini<n; ini+=L_ini){
      // initial superimposed fragment
      int end=ini+L_ini; if(end>n)break;
      for(k=0; k<n; k++){
	if((k>=ini)&&(k<end)&&(mprot[k]>1)){w[k]=1;}else{w[k]=0;}
      }

      float TM_max_tmp=0, TM=0; int nfail=0;
      for(int iter=0; iter<n; iter++){
	int align=0; for(k=0; k<n; k++)if(w[k])align++;
	if(align<nali_min)break;

	if(Superimpose_mult(w, xca, len, mprot, i_ref, msa, nprot, n)<0){
	  break;
	}
	for(i=0; i<nprot; i++){ // Check if nan, in case restore
	  for(j=0; j<ncart[i]; j++)
	    if(isnan(xca[i][j])){
	      Copy_coord(xca[i], prot_p[i]->xca_rot, ncart[i]); break;
	    }
	}

	TM=Compute_TM_mult(tm_k, NULL, xca, msa, nprot, n, d02, len);
	//printf("TM: nali=%d %.1f nfail=%d\n", align, TM/lmin, nfail);

	if(TM>TM_max_tmp){
	  TM_max_tmp=TM; nfail=0;
	  if(TM>TM_max){TM_max=TM; for(i=0; i<n; i++)w_opt[i]=w[i];}
	}else{
	  nfail++; if(nfail > nfail_max){break;}
	}
	
	// Superimpose only if tm_k > thr
	float tm_ave=TM/n;
	for(i=0; i<n; i++){
	  float thr=(iter*tm_thr[i]+tm_ave)/(1+iter);
	  if(mprot[i]>1 && tm_k[i]>thr){w[i]=1;}else{w[i]=0;}
	}
      } // end rounds

    } // end
  } // end

  for(i=0; i<nprot; i++){ // Check if nan, in case restore
    for(j=0; j<ncart[i]; j++)
      if(isnan(xca[i][j])){
	Copy_coord(xca[i], prot_p[i]->xca_rot, ncart[i]); break;
      }
  }

  double TM;
  if(Superimpose_mult(w_opt, xca, len, mprot, i_ref, msa, nprot, n)<0){
    printf("ERROR, multiple superimposition did not work\n");
    TM=-1;
  }else{
    TM=Compute_TM_mult(tm_k, TM_all, xca, msa, nprot, n, d02, len);
    for(i=0; i<nprot; i++)Copy_coord(prot_p[i]->xca_rot, xca[i], ncart[i]);
  }
  for(i=0; i<nprot; i++)free(xca[i]);

  return(TM/lmin);
}

float Compute_TM_mult(double *tm_k, float **TM_all, float **xca,
		      int **msa, int nprot, int n, float d02, int *len)
{
  int k; for(k=0; k<n; k++){tm_k[k]=0;}
  for(int i=0; i<nprot; i++){
    int *msa_i=msa[i]; float *xca_i=xca[i];
    for(int j=0; j<i; j++){
      int *msa_j=msa[j]; float *xca_j=xca[j];
      double TM_ij=0;
      for(k=0; k<n; k++){
	int ki=msa_i[k]; if(ki<0)continue;
	int kj=msa_j[k]; if(kj<0)continue;
	float dd=dist2(xca_i+3*ki, xca_j+3*kj);
	float tm=1./(1+dd/d02);
	tm_k[k]+=tm; if(TM_all)TM_ij+=tm;
      }
      if(TM_all){
	int ltar;
	if(len[i]<len[j]){ltar=len[i];}else{ltar=len[j];}
	TM_all[i][j]=TM_ij/ltar;
      }
    }
  }
  double TM=0;
  for(k=0; k<n; k++){TM+=tm_k[k];}
  if(isnan(TM)){
    printf("ERROR, TM is nan d0=%.2g\ntm_k= ", d02);
    for(k=0; k<n; k++)printf(" %.3g", tm_k[k]);
    printf("\n"); //exit(8);
  }
  return(TM);
}

int Superimpose_mult(float *w, float **xca, int *len, int *mprot,
		      int i_ref, int **msa, int nprot, int n)
{
  int nmin2=nali_min/2, ret=0;
  int nn=0, iali[n], k, j;
  for(k=0; k<n; k++){
    if(w[k]==0)continue;
    iali[nn]=k; nn++;
  }
  if(nn<nali_min){
    printf("WARNING, too few aligned residues (%d) in Superimpose_mult\n",
	   nn); return(-1);
  }

  // Initialize average conformation
  int sum[nn]; float xsum[3*nn], *xa=xsum;
  for(k=0; k<nn; k++){
    sum[k]=0; for(int j=0; j<3; j++){*xa=0; xa++;}
  }
  Add_coord(xsum, sum, nn, xca[i_ref], msa[i_ref], iali);
  
  // align to average conformation
  float *rot[3], shift[3];
  for(int i=0; i<3; i++)rot[i]=malloc(3*sizeof(float));

  float xca1[3*nn], xca2[3*nn];
  for(int iter=0; iter<N_SUP_AVE; iter++){
    int aligned=1;
    for(int i=0; i<nprot; i++){
      float *xca_i=xca[i]; int *msa_i=msa[i];
      if(iter){Remove_coord(xsum, sum, nn, xca_i, msa_i, iali);}
      else{if(i==i_ref)continue;}
      int nali=0; float *x1=xca1, *x2=xca2, ww[nn];
      for(k=0; k<nn; k++){
	int kk=iali[k];
	if(kk>n)printf("WARNING prot %d k=%d kk= %d\n", i, k, kk);
	if(sum[k]==0 || msa_i[kk]<0)continue;
	float *xs=xsum+3*k, *xi=xca_i+3*msa_i[kk];
	for(j=0; j<3; j++){
	  *x1=*xs; if(sum[k]){*x1/=sum[k];} xs++; x1++;
	  *x2=*xi; xi++; x2++;
	}
	ww[nali]=sum[k]; nali++;
      }
      if(nali>nmin2){
	rmsd_mclachlan_f(rot, shift, xca1, xca2, ww, nali);
	Rot_shift(xca_i, len[i], rot, shift);
	Add_coord(xsum, sum, nn, xca_i, msa_i, iali);
	aligned++;
      }else{
	//printf("WARNING, too few aligned residues (%d) for str %d"
	//     " iter %d nn=%d\n", nali, i, iter,nn);
      }
    }
    if(aligned<nprot){
      ret=-1;
      printf("WARNING, %d proteins could not be aligned iter=%d nali=%d\n",
	     nprot-aligned, iter, nn);
    }else{
      ret=0; //if(iter==1)goto empty_rot;
    }
  } // end iter
 empty_rot:
  for(int i=0; i<3; i++)free(rot[i]);
  return(ret);
}

void Rot_shift(float *xca_i, int len, float **rot, float *shift)
{
  float *x=xca_i;
  for(int k=0; k<len; k++){
    if(rot)Rotate_vector(x, rot);
    if(shift){for(int j=0; j<3; j++){*x+=shift[j]; x++;}}
    else{x+=3;}
  }
}

void Add_coord(float *xsum, int *sum, int nn, float *xca_i, int *msa_i,
	       int *iali)
{
  for(int k=0; k<nn; k++){
    int k1=msa_i[iali[k]];
    if(k1>=0){
      float *xca_k=xca_i+3*k1, *xs_k=xsum+3*k;
      for(int j=0; j<3; j++){*xs_k+=*xca_k; xs_k++; xca_k++;}
      sum[k]++;
    }
  }
}

void Remove_coord(float *xsum, int *sum, int nn, float *xca_i, int *msa_i,
		  int *iali)
{
  for(int k=0; k<nn; k++){
    int k1=msa_i[iali[k]];
    if(k1>=0){
      float *xca_k=xca_i+3*k1, *xs_k=xsum+3*k;
      for(int j=0; j<3; j++){*xs_k-=*xca_k; xs_k++; xca_k++;}
      sum[k]--;
    }
  }
}


float TM_fast(float **d2_out, float d02, int *ali, int ltar, float *d2min1,
	      struct protein *prot1, struct protein *prot2)
{
  float *xca1_store=prot1->xca_rot; int n1=prot1->len;
  float *xca2_store=prot2->xca_rot; int n2=prot2->len;

  float TM_max=0, w[n1], w_TM[n1];
  float xca1[3*n1], xca2[3*n2];
  int n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
  int i, k=0, align=0;
  for(i=0; i<n1; i++){
    if(ali[i]<0)continue;
    if(d2min1[i]<d02){w[k]=1; align++;}else{w[k]=0;} k++;
  }
  for(int round=0; round<n; round++){
    if(align<4)break;
    rmsd_mclachlan_f(NULL, NULL, xca1, xca2, w, n);
    double TM=0, d2_max=0; int k=0, imax=-1;
    for(i=0; i<n; i++){
      float d2=dist2(xca1+k, xca2+k); k+=3;
      TM+=1./(1+d2/d02);
      if(w[i]&&(d2>d2_max)){d2_max=d2; imax=i;}
    }
    if(TM>TM_max){TM_max=TM; for(i=0; i<n; i++)w_TM[i]=w[i];}
    else {break;}
    // Remove most divergent pair
    if(imax>=0){w[imax]=0; align--;}
  }
   // Superimpose all residues, including not aligned
  float shift[3], **rot=Allocate_mat2_f(3, 3); int j;
  n=Set_coord_ali(xca1, xca1_store, xca2, xca2_store, ali, n1);
  rmsd_mclachlan_f(rot, shift, xca1, xca2, w_TM, n);
  Rot_shift(xca2_store, n2, rot, shift);
  Empty_matrix_f(rot, 3);

 if(d2_out){
    All_distances(d2_out, xca1_store, n1, xca2_store, n2);
  }
  return(TM_max/ltar);
}

void Examine_neighbors(float **d2, int *ali, float d02,
		       int *shift, int *id_3D,
		       int *neigh_ali, int *neigh_noali,
		       int *neigh_noali_aaid,
		       char *seq1, int n1, char *seq2, int n2)
{
  // fraction of d0 to identify superimposed residues
  float d_thr=0.5, thr=d_thr*d02; int i;

  // Identify aligned and well superimposed residues below threshold
  if(id_3D){
    for(i=0; i<n1; i++){
      if(ali[i]<0){id_3D[i]=0; continue;}
      else{
	if(d2[i][ali[i]]<thr){id_3D[i]=1;}else{id_3D[i]=0;}
      }
    }
  }

  // Find the closest atom of str.2 for each atom of str.1 and vice-versa
  int neighb1[n1]; float d2min1[n1], score1[n1];
  Closest_neighbor(neighb1, d2min1, score1, d02, d2, n1, n2, 0);
  int neighb2[n2]; float d2min2[n2], score2[n2];
  Closest_neighbor(neighb2, d2min2, score2, d02, d2, n2, n1, 1);

  // Find neighbors that are not aligned
  if(neigh_noali && neigh_ali){
    int ali2[n2]; Invert_ali(ali2, n2, ali, n1);
    for(i=0; i<n1; i++){
      int j=neighb1[i];
      if((j>=0)&&(neighb2[j]==i)&&(d2min1[i]<thr)){
	if(ali[i]==j){neigh_ali[i]=1;}else{neigh_noali[i]=1;}
      }else{
	neigh_ali[i]=0;
	neigh_noali[i]=0;
      }

      if(id_3D[i] && neigh_ali[i]==0 && j>=0){
	int s1=abs(i-neighb2[j]), s2=abs(j-ali[i]);
	if(s1>s2){shift[i]=s1;}else{shift[i]=s2;}
      }else{
	shift[i]=0;
      }

      if(neigh_noali[i]){
	// Residue i is closest to not aligned residue
	// is ali identical?
	if(seq1 && ((ali[i] >=0 && seq1[i]==seq2[ali[i]])||
		    (j>=0 && ali2[j]>=0 && seq2[j]==seq1[ali2[j]]))){
	  neigh_noali_aaid[i]=1;	  
	}else{
	  neigh_noali_aaid[i]=0;
	}
      }
    }
  }

}

void Align_TM(int *ali_new, float *d2min1, float **d2, float d02, int *ali,
	      int n1, int n2)
{
  // Change alignment, optimize the TM score of new alignment
  float tthr=1./(1+0.75);

  // Find the closest atom of str.2 for each atom of str.1 and vice-versa
  int neighb1[n1]; float score1[n1]; //d2min1[n1], 
  Closest_neighbor(neighb1, d2min1, score1, d02, d2, n1, n2, 0);
  int neighb2[n2]; float d2min2[n2], score2[n2];
  Closest_neighbor(neighb2, d2min2, score2, d02, d2, n2, n1, 1);
  Align_neighbors(ali_new, ali, tthr, neighb1,score1,n1,neighb2,score2,n2);
}


int Set_weights(int *w1, int *w2, float *w_ali, int *ali, int n1, int n2)
{
  int i, k=0;
  for(i=0; i<n1; i++)w1[i]=0;
  for(i=0; i<n2; i++)w2[i]=0;
  for(i=0; i<n1; i++){
    if(ali[i]<0)continue;
    if(w_ali[k]){w1[i]=1; w2[ali[i]]=1;} k++;
  }
  return(0);
}

void All_distances(float **d2, float *x1, int n1, float *x2, int n2)
{
  float *x1i=x1;
  for(int i=0; i<n1; i++){
    float *x2j=x2;
    for(int j=0; j<n2; j++){
      d2[i][j]=dist2(x1i, x2j); x2j+=3;
    }
    x1i+=3;
  }
}

void Closest_neighbor(int *neighb1, float *dmin1, float *score, float d0,
		      float **d2, int n1, int n2, int order)
{
  for(int i=0; i<n1; i++){
    float d_min=1000, dd; int j_min=-1;
    for(int j=0; j<n2; j++){
      if(order==0){dd=d2[i][j];}else{dd=d2[j][i];}
      if(dd<d_min){d_min=dd; j_min=j;}
    }
    neighb1[i]=j_min; dmin1[i]=d_min; score[i]=1./(1+d_min/d0);
  }
}


void Rotate_vector(float *x, float **rot){
  float v[3];
  for(int i=0; i<3; i++){
    double w=0; for(int j=0; j<3; j++)w+=rot[i][j]*x[j]; v[i]=w;
  }
  for(int i=0; i<3; i++)x[i]=v[i];
}

int Set_coord_cm(float *x, double *cm, float *x_store, int *w, int n){
  int i, j, m=0;
  for(i=0; i<3; i++)cm[i]=0;
  float *xj=x, *xs=x_store;
  for(i=0; i<n; i++){
    if(w[i])m++;
    for(j=0; j<3; j++){*xj=*xs; if(w[i])cm[j]+=*xs; xs++; xj++;}
  }
  for(i=0; i<3; i++)cm[i]/=m;
  xj=x;
  for(i=0; i<n; i++){
    for(j=0; j<3; j++){*xj-=cm[j]; xj++;}
  }
  return(0);
}

float dist2(float *x1, float *x2){
  float dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
  return(dx*dx+dy*dy+dz*dz);
}

int Set_coord_ali(float *xca1, float *xca1_store,
		  float *xca2, float *xca2_store, int *ali, int n1)
{
  // Set coordinates of protein 2 aligned to protein 1
  int n=0; 
  float *x1=xca1, *x2=xca2;
  for(int i1=0; i1<n1; i1++){
    int i2=ali[i1]; if(i2<0)continue;
    n++;
    float *xs1=xca1_store+3*i1, *xs2=xca2_store+3*i2;
    for(int j=0; j<3; j++){
      *x1=*xs1; x1++; xs1++;
      *x2=*xs2; x2++; xs2++;
    }
  }
  return(n);
}

void Copy_coord(float *xca, float *xca_store, int ncart)
{
  float *xj=xca, *xs=xca_store;
  for(int i=0; i<ncart; i++){*xj=*xs; xj++; xs++;}
}


void Test_Rot(float **d2, int iprot, int jprot,
	      float **Rot, float *Shift,
	      float *xca1_store, int len1,
	      float *xca2_store, int len2)
{
  // Rotate structure 2
  float xca1[3*len1], xca2[3*len2], *d2_test[len1]; int i;
  Copy_coord(xca1, xca1_store, 3*len1);
  Copy_coord(xca2, xca2_store, 3*len2);
  Rot_shift(xca2, len2, Rot, Shift);
  for(i=0; i<len1; i++)d2_test[i]=malloc(len2*sizeof(float));
  All_distances(d2_test, xca1, len1, xca2, len2);
 
  int err=0; float *x=xca1;
  for(i=0; i<len1; i++){
    if(isnan(*x) || isnan(*(x+1)) || isnan(*(x+2))){
      float *xs=xca1_store+3*i;
      printf("ERROR in coord1 %d str. %d-%d L=%d "
	     "r= %.2g %.2g %.2g r_store= %.2g %.2g %.2g\n",
	     i, iprot, jprot, len1, *x, *(x+1), *(x+2),
	     *(xs), *(xs+1), *(xs+2));
      err++;
    }
    x+=3;
  }
  err=0; x=xca2;
  for(i=0; i<len2; i++){
    if(isnan(*x) || isnan(*(x+1)) || isnan(*(x+2))){
      float *xs=xca2_store+3*i;
      printf("ERROR in coord2 %d str. %d-%d L=%d "
	     "r= %.2g %.2g %.2g r_store= %.2g %.2g %.2g\n",
	     i, jprot, iprot, len2, *x, *(x+1), *(x+2),
	     *(xs), *(xs+1), *(xs+2));
      err++;
    }
    x+=3;
  }
  if(err){
    printf("Shift: %.2g %.2g %.2g\n", Shift[0], Shift[1], Shift[2]);
    printf("Rot: %.2g %.2g %.2g %.2g %.2g %.2g %.2g %.2g %.2g\n", 
	   Rot[0][0],Rot[0][1],Rot[0][2],
	   Rot[1][0],Rot[1][1],Rot[1][2],
	   Rot[2][0],Rot[2][1],Rot[2][2]);
  }


  double diff_d2=0;
  for(i=0; i<len1; i++){
    for(int j=0; j<len2; j++)diff_d2+=d2[i][j]-d2_test[i][j];
  }
  diff_d2=sqrt(diff_d2/(len1*len2));
  if(diff_d2 > 0.1){
    printf("ERROR, iprot= %d jprot= %d sqrt(<(d^2-d_test^2)>)= %.2g\n",
	   iprot, jprot, diff_d2); exit(8);
  }

  if(err)exit(8);
  for(i=0; i<len1; i++)free(d2_test[i]);
}

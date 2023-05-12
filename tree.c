#include "tree.h"
#include <stdio.h>
#include <stdlib.h>

int Write_cluster(int **elements, int *clus, int n_clus,
		  int *cluster, int n, float **div, int ini);


struct treenode *Neighbor_Joining(float **div, int n){
  int N_nodes=2*n-2, i;
  struct treenode *node=malloc(N_nodes*sizeof(struct treenode));
  for(i=0; i<n; i++){
    node[i].offs1=NULL;
    node[i].offs2=NULL;
    node[i].parent=NULL;
  }
  return(node);
}

int ***Build_outgroups(int ***N_out, struct treenode *node, int n)
{
  int ***outgroup=malloc(n*sizeof(int **));
  *N_out=malloc(n*sizeof(int *));
  int out[n], i, j, k;
  for(i=0; i<n; i++){
    (*N_out)[i]=malloc(i*sizeof(int));
    outgroup[i]=malloc(i*sizeof(int *));
    for(j=0; j<i; j++){
      int m=0; // number of outgroups
      for(k=0; k<n; k++){
	if((k==i)||(k==j))continue;
	if(0){ // Condition for outgroups
	  out[m]=k; m++;
	}
      }
      (*N_out)[i][j]=m;
      int *out2=malloc(m*sizeof(int));
      for(k=0; k<m; k++)out2[k]=out[k];
      outgroup[i][j]=out2;
    }
  }
  return(outgroup);
}

int ***Outgroups_NJ(int ***N_out, float **div, int n)
{
  int ***outgroup=malloc(n*sizeof(int **));
  *N_out=malloc(n*sizeof(int *));
  int out_tmp[n], i, j, k;

  // Average distance
  float ave_div[n]; int n2=n-2;
  for(i=0; i<n; i++){
    double q=0;
    for(j=0; j<n; j++)if(i!=j)q+=div[i][j];
    ave_div[i]=q/n2;
  }

  for(i=0; i<n; i++){
    (*N_out)[i]=malloc((i-1)*sizeof(int));
    outgroup[i]=malloc((i-1)*sizeof(int *));
    for(j=0; j<i; j++){
      int m=0; // number of outgroups
      float q_ij=div[i][j]-ave_div[j], q_ji=div[i][j]-ave_div[i];
      for(k=0; k<n; k++){
	if((k==i)||(k==j))continue;
	if((q_ij<(div[i][k]-ave_div[k]))&&
	   (q_ji<(div[j][k]-ave_div[k]))){
	  out_tmp[m]=k; m++;
	}
      }
      (*N_out)[i][j]=m;
      int *outg=malloc(m*sizeof(int));
      outgroup[i][j]=outg;
      for(k=0; k<m; k++)outg[k]=out_tmp[k];
    }
  }
  return(outgroup);
}

int Single_linkage(int ***elements, int **n_ele, int **clus,
		   float **div, int n, float thr)
{
  // Make single linkage from the point of view of elements
  int nclus=n,  cluster[n], i;
  for(i=0; i<n; i++)cluster[i]=i;
  int joined=1, round=0, c_old, c_new;
  while(joined){
    round++; printf("Joining clusters, round %d\n", round);
    joined=0;
    for(i=0; i<n; i++){
      for(int j=0; j<i; j++){
	if((cluster[i]!=cluster[j])&&(div[i][j]<thr)){
	  joined=1; nclus--; // Join clusters i and j
	  if(cluster[i]<cluster[j]){c_new=cluster[i]; c_old=cluster[j];}
	  else{c_new=cluster[j]; c_old=cluster[i];}
	  for(int k=0; k<n; k++)
	    if(cluster[k]==c_old)cluster[k]=c_new;
	}
      }
    }
  }
  //  returned: clus[i], nclus

  // Store properties of the clusters
  *n_ele=malloc(nclus*sizeof(int));      // Number of elements in a cluster
  *elements=malloc(nclus*sizeof(int *)); // List of cluster elements
  *clus=malloc(n*sizeof(int));           // Cluster index of element i
  for(i=0; i<n; i++)(*clus)[i]=-1;
  int n_clus=0;
  for(i=0; i<n; i++){
    if(cluster[i]>=0){
      (*n_ele)[n_clus]=
	Write_cluster(*elements+n_clus, *clus, n_clus, cluster, n, div, i);
      n_clus++;
    }
  }
  if(n_clus != nclus){
    printf("ERROR in Single_likage, wrong number of clusters %d exp. %d\n",
	   n_clus, nclus); exit(8);
  }
  for(i=0; i<n; i++){
    if(((*clus)[i]<0)||((*clus)[i]>=n_clus)){
      printf("ERROR in cluster index\n"); exit(8);
    }
  }

  return(n_clus);
}

int Write_cluster(int **elements, int *clus, int n_clus,
		  int *cluster, int n, float **div, int ini)
{
  // Write cluster
  int new=cluster[ini]; // label of cluster
  int m=0; // index of element in cluster
  int c_ele[n], j; float d_sum[n];
  for(j=ini; j<n; j++){
    if(cluster[j]==new){
      clus[j]=n_clus;
      cluster[j]=-1;
      c_ele[m]=j;
      d_sum[m]=0;
      for(int k=0; k<m; k++){
	int l=c_ele[k];
	d_sum[m]+=div[j][l];
	d_sum[k]+=div[j][l];
      }
      m++;
    }
  }
  (*elements)=malloc(m*sizeof(int));
  float d_min=d_sum[0];
  for(j=0; j<m; j++){
    (*elements)[j]=c_ele[j];
    if(d_sum[j]<d_min){d_min=d_sum[j];}
  }
  return(m);
}

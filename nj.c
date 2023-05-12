// Neighbor joining

struct node{
  struct node *parent, *son1, *son2;
  float branch;
  int nson, *sons;
  char *name;
  int index;
};

void Copy_sons(struct node *mnode, struct node *node1, struct node *node2);
void Copy_name(struct node *mnode, struct node *node1, struct node *node2);

#include "protein.h"
#include "nj_align.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int NJ_align(int **msa, int *N_ali_max,
	     int ***Ali_pair, float ****Rot_pair, float ***Shift_pair,
	     struct protein **prot_p, char *name_out, 
	     float **div_in, int n, int AV_LINK)
{
  printf("Building tree of %d proteins with ", n);
  if(AV_LINK)printf("Average linkage");
  else{printf("Neighbor Joining");}
  printf(" method\n");

  int L_msa=0; 
  int n_nodes=2*n-1, a, b;
  double Dsum[n_nodes];
  float *div[n_nodes], *Q[n_nodes];
  for(a=0; a<n_nodes; a++){
    div[a]=malloc(n_nodes*sizeof(float));
  }
  // Initialize matrices for neighbor joining
  for(a=0; a<n; a++)Dsum[a]=0;
  for(a=0; a<n; a++)
    for(b=a+1; b<n; b++){
      float d=div_in[a][b];
      div[a][b]=d; div[b][a]=d;
      Dsum[a]+=d; Dsum[b]+=d;
    }
  int nout=n-2; // number of outgroups
  if(AV_LINK==0){
    for(a=0; a<n_nodes; a++)Q[a]=malloc(n_nodes*sizeof(float));
    for(a=0; a<n; a++){
      for(b=a+1; b<n; b++){
	Q[a][b]=nout*div[a][b]-Dsum[a]-Dsum[b];
	Q[b][a]=Q[a][b];
      }
    }
  }

  // Initialize tree
  int label[n]; for(a=0; a<n; a++)label[a]=a;
  struct node nodes[n_nodes], *node=nodes;
  for(a=0; a<n_nodes; a++){
    node->index=a;
    node->son1=NULL; node->son2=NULL; node->parent=NULL;
    if(a<n){
      node->nson=1;
      node->sons=malloc(1*sizeof(int)); node->sons[0]=a;
      node->name=(char *)malloc(sizeof(prot_p[a]->name));
      strcpy(node->name, prot_p[a]->name);
    }else{
      node->nson=0; node->sons=NULL; node->name=NULL;
    }
    node->branch=0;
    node++;
  }

  // Initialize MSA
  int len_msa[n], i;
  if(msa){
    for(a=0; a<n; a++){
      len_msa[a]=prot_p[a]->len;
      for(i=0; i< len_msa[a]; i++)msa[a][i]=i;
    }
  }

  float **score; if(AV_LINK){score=div;}else{score=Q;}
  struct node *mnode; int mum=n; // new parent node
  for(int nclus=n; nclus>1; nclus--){
    // Look for minimum of Q
    int amin=-1, lamin, bmin, lbmin;
    float Qmin=10000000, d;
    for(a=0; a<nclus; a++){
      int la=label[a];
      for(b=a+1; b<nclus; b++){
	int lb=label[b];
	if(score[la][lb]<Qmin || amin<0){
	  Qmin=score[la][lb]; amin=a; lamin=la; bmin=b; lbmin=lb;
	}
      }
    }

    // Join lamin and lbmin
    printf("Joining nodes %d %d into node %d -> %d clusters\n",
	   lamin, lbmin, mum, nclus-1);
    mnode=nodes+mum;
    nodes[lamin].parent=mnode;
    nodes[lbmin].parent=mnode;
    mnode->son1=nodes+lamin;
    mnode->son2=nodes+lbmin;
    Copy_sons(mnode, mnode->son1, mnode->son2);
    // Distance from mother node
    if(nout){
      d=0.5*(div[lamin][lbmin]+(Dsum[lamin]-Dsum[lbmin])/nout);
    }else{
      d=0.5*div[lamin][lbmin];
    }
    nodes[lamin].branch=d;
    nodes[lbmin].branch=div[lamin][lbmin]-d;
    Copy_name(mnode, mnode->son1, mnode->son2);
    printf("%s\n", mnode->name);
    // Align amin and bmin
    if(msa){
      L_msa=Align_PC_NW_mult(msa, len_msa, N_ali_max, n, prot_p,
			     Ali_pair, Rot_pair, Shift_pair,
			     mnode->son1->sons, mnode->son1->nson,
			     mnode->son2->sons, mnode->son2->nson);
    }

    // Outliers
    Dsum[mum]=0;
    for(a=0; a<nclus; a++){
      if(a==amin || a==bmin)continue;
      int la=label[a];
      float d=div[lamin][la]+div[lbmin][la];
      if(AV_LINK){
	div[la][mum]=
	  (nodes[lamin].nson*div[lamin][la]+nodes[lbmin].nson*div[lbmin][la])
	  /(nodes[lamin].nson+nodes[lbmin].nson);
      }else{
	div[la][mum]=0.5*(d-div[lamin][lbmin]);
      }
      div[mum][la]=div[la][mum];
      Dsum[la]+=(div[la][mum]-d);
      Dsum[mum]+=div[la][mum];
    }

    // update the Q matrix
    nout--;
    if(AV_LINK==0){
      for(a=0; a<nclus; a++){
	if(a==amin || a==bmin)continue;
	int la=label[a];	
	Q[mum][la]=nout*div[mum][la]-Dsum[mum]-Dsum[la];
	Q[la][mum]=Q[mum][la];
	for(b=a+1; b<nclus; b++){
	  int lb=label[b];
	  if(b==amin || b==bmin)continue;
	  Q[la][lb]=nout*div[la][lb]-Dsum[la]-Dsum[lb];
	  Q[lb][la]=Q[la][lb];
	}
      }
    }

    // Update labels
    label[amin]=mum;
    for(b=bmin; b<(nclus-1); b++)label[b]=label[b+1];
    mum++;

  } // end clustering

  // Write the tree in Newick format
  char nametree[100]; sprintf(nametree, "%s.tree", name_out);
  printf("Writing tree in Newick format in %s\n", nametree);
  FILE *file_out=fopen(nametree, "w");
  fprintf(file_out, "%s\n", mnode->name);
  fclose(file_out);

  for(a=0; a<n_nodes; a++){
    free(div[a]); if(AV_LINK==0)free(Q[a]);
  }
  node=nodes;
  for(a=0; a<n_nodes; a++){
    free(node->sons); free(node->name); node++;
  }
  return(L_msa);
}

void Copy_sons(struct node *mnode, struct node *node1, struct node *node2)
{
  int k, i;
  mnode->nson=node1->nson+node2->nson;
  mnode->sons=malloc(mnode->nson*sizeof(int));
  for(i=0; i<node1->nson; i++){mnode->sons[i]=node1->sons[i];} k=i;
  for(i=0; i<node2->nson; i++){mnode->sons[k]=node2->sons[i]; k++;}
}

void Copy_name(struct node *mnode, struct node *node1, struct node *node2)
{
  mnode->name=
    malloc((strlen(node1->name)+strlen(node2->name)+30)*sizeof(char));
  sprintf(mnode->name, "(%s:%.3g,%s:%.3g)",
	  node1->name, node1->branch, node2->name, node2->branch);
}

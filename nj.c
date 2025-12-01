// Neighbor joining

#include <math.h>

struct node{
  int index;
  struct node *parent, *son1, *son2;
  int nson, *sons;
  float branch, branch_clock;
  char *name;
  float branch_seq, branch_seq_clock; // based on sequence
  char *name_seq;
};

void Copy_sons(struct node *mnode, struct node *node1, struct node *node2);
void Copy_name(struct node *mnode, struct node *node1, struct node *node2);
char *Copy_name_3(struct node *node1, struct node *node2, struct node *node3);
void Copy_name_seq(struct node *mnode, struct node *node1, struct node *node2);
char *Copy_name_3_seq(struct node *node1,struct node *node2,struct node *node3);


#include "protein.h"
#include "nj_align.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int NJ_align(int **msa, int *N_ali_max,
	     int ***Ali_pair, float ****Rot_pair, float ***Shift_pair,
	     struct protein **prot_p, char *name_out, 
	     float **div_in, float **div_seq_in,
	     int n, int AV_LINK)
{
  int NJ=0, NJ_BR=1;
  printf("Building tree of %d proteins with ", n);
  if(AV_LINK){printf("Average linkage"); NJ=0;}
  else{printf("Neighbor Joining"); NJ=1;}
  printf(" method\n");
  if(NJ){NJ_BR=1;} // determine  branch lengths with mean square method
  if(NJ_BR){
    printf("The branch lengths are determined by fitting the distances"
	   " (no molecular clock assumed)\n");
  }

  int L_msa=0; 
  int n_nodes=2*n-1, a, b;
  double Dsum[n_nodes];
  float *div[n_nodes], *Q[n_nodes];
  for(a=0; a<n_nodes; a++){
    div[a]=malloc(n_nodes*sizeof(float));
  }

  // Initialize matrices for neighbor joining
  for(a=0; a<n; a++)Dsum[a]=0;
  for(a=0; a<n; a++){
    for(b=a+1; b<n; b++){
      float d=div_in[a][b];
      div[a][b]=d; div[b][a]=d;
      if(NJ_BR){Dsum[a]+=d; Dsum[b]+=d;}
    }
  }

  // Sequence divergence
  float **div_seq=NULL; double Dsum_seq[n_nodes];
  if(div_seq_in){
    div_seq=malloc(n_nodes*sizeof(float *));
    for(a=0; a<n_nodes; a++){
      div_seq[a]=malloc(n_nodes*sizeof(float));
    }
    for(a=0; a<n; a++)Dsum_seq[a]=0;
    for(a=0; a<n; a++){
      div_seq[a][a]=0;
      for(b=a+1; b<n; b++){
	float d=div_seq_in[a][b];
	div_seq[a][b]=d; div_seq[b][a]=d;
	if(NJ_BR){Dsum_seq[a]+=d; Dsum_seq[b]+=d;}
      }
    }
  }

  int nout=n-2; // number of outgroups
  if(NJ){
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
      node->name=(char *)malloc(sizeof(prot_p[a]->domname));
      strcpy(node->name, prot_p[a]->domname);
      if(div_seq){
	node->name_seq=(char *)malloc(sizeof(prot_p[a]->domname));
	strcpy(node->name_seq, prot_p[a]->domname);
      }else{
	node->name_seq=NULL;
      }
    }else{
      node->nson=0; node->sons=NULL;
      node->name=NULL; node->name_seq=NULL;
    }
    node->branch=0;
    node->branch_seq=0;
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

  float **score; if(NJ){score=Q;}else{score=div;}
  struct node *mnode;  int mum=n; // new parent node
  char *umnode=NULL, *umnode_seq=NULL;  // unrooted
  for(int nclus=n; nclus>1; nclus--){
    // Look for minimum of Q
    int amin=-1, lamin, bmin, lbmin;
    float Smin=10000000;
    for(a=0; a<nclus; a++){
      int la=label[a];
      for(b=a+1; b<nclus; b++){
	int lb=label[b];
	if(score[la][lb]<Smin || amin<0){
	  Smin=score[la][lb]; amin=a; lamin=la; bmin=b; lbmin=lb;
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
    float d=0.5*div[lamin][lbmin], dd=0;
    if(NJ_BR && nout){
      dd=0.5*(Dsum[lamin]-Dsum[lbmin])/nout;
      if(dd>d || dd<-d){dd=0.5*d;}
    }
    nodes[lamin].branch=d-dd;
    nodes[lbmin].branch=d+dd;
    Copy_name(mnode, mnode->son1, mnode->son2);
    if(div_seq){
      float d=0.5*div_seq[lamin][lbmin], dd=0;
      if(NJ_BR && nout){
	dd=0.5*(Dsum_seq[lamin]-Dsum_seq[lbmin])/nout;
	if(dd>d || dd<-d){dd=0.5*d;}
      }
      nodes[lamin].branch_seq=d-dd;
      nodes[lbmin].branch_seq=d+dd;
      Copy_name_seq(mnode, mnode->son1, mnode->son2);
    }
    //printf("%s\n", mnode->name); // change: commented out
    // Align amin and bmin
    if(msa){
      L_msa=Align_PC_NW_mult(msa, len_msa, N_ali_max, n, prot_p,
			     Ali_pair, Rot_pair, Shift_pair,
			     mnode->son1->sons, mnode->son1->nson,
			     mnode->son2->sons, mnode->son2->nson);
    }

    // Distance between mother node and outgroups
    Dsum[mum]=0;
    Dsum_seq[mum]=0;
    for(int o=0; o<nclus; o++){
      if(o==amin || o==bmin)continue;
      int lo=label[o];
      if(NJ_BR){
	float d=div[lamin][lo]+div[lbmin][lo];
	div[lo][mum]=0.5*(d-div[lamin][lbmin]);
	Dsum[lo]+=(div[lo][mum]-d);
	Dsum[mum]+=div[lo][mum];
      }else{
	div[lo][mum]=
	  (nodes[lamin].nson*div[lamin][lo]+nodes[lbmin].nson*div[lbmin][lo])
	  /(nodes[lamin].nson+nodes[lbmin].nson);
      }
      div[mum][lo]=div[lo][mum];

      if(div_seq){
	if(NJ_BR){
	  float d=div_seq[lamin][lo]+div_seq[lbmin][lo];
	  div_seq[lo][mum]=0.5*(d-div_seq[lamin][lbmin]);
	  Dsum_seq[lo]+=(div_seq[lo][mum]-d);
	  Dsum_seq[mum]+=div_seq[lo][mum];
	}else{
	  div_seq[lo][mum]=
	    (nodes[lamin].nson*div_seq[lamin][lo]+
	     nodes[lbmin].nson*div_seq[lbmin][lo])
	    /(nodes[lamin].nson+nodes[lbmin].nson);
	}
	div_seq[mum][lo]=div_seq[lo][mum];
      }
    }
    
    if(nclus==3){
      // unrooted tree
      int lc=-1;
      for(int o=0; o<nclus; o++){
	int lo=label[o];
	if(lo != lamin && lo != lbmin){lc=lo; break;}
      }
      if(lc>=0){
	nodes[lc].branch=div[mum][lc];
	umnode=Copy_name_3(mnode->son1, mnode->son2, nodes+lc);
	if(div_seq){
	  nodes[lc].branch_seq=div_seq[mum][lc];
	  umnode_seq=Copy_name_3_seq(mnode->son1, mnode->son2, nodes+lc);
	}
      } 
    }

    // update the Q matrix
    nout--;
    if(NJ){
      for(int o=0; o<nclus; o++){
	if(a==amin || a==bmin)continue;
	int lo=label[o];	
	Q[mum][lo]=nout*div[mum][lo]-Dsum[mum]-Dsum[lo];
	Q[lo][mum]=Q[mum][lo];
	for(int p=o+1; p<nclus; p++){
	  int lp=label[p];
	  if(p==amin || p==bmin)continue;
	  Q[lo][lp]=nout*div[lo][lp]-Dsum[lo]-Dsum[lp];
	  Q[lp][lo]=Q[lo][lp];
	}
      }
    }

    // Update labels
    label[amin]=mum;
    for(b=bmin; b<(nclus-1); b++)label[b]=label[b+1];
    mum++;

  } // end clustering

  // Write the tree in Newick format
  char nametree[200]; FILE *file_out;
  if(0){
    if(NJ){sprintf(nametree, "%s.rooted.nj.tree", name_out);}
    else{sprintf(nametree, "%s.rooted.al.tree", name_out);}
    printf("Writing rooted tree in Newick format in %s\n", nametree);
    file_out=fopen(nametree, "w");
    fprintf(file_out, "%s\n", mnode->name);
    fclose(file_out);
    if(div_seq){
      if(NJ){sprintf(nametree, "%s.rooted.nj.seqtree", name_out);}
      else{sprintf(nametree, "%s.rooted.al.seqtree", name_out);}
      printf("Writing rooted seq tree in Newick format in %s\n", nametree);
      file_out=fopen(nametree, "w");
      fprintf(file_out, "%s\n", mnode->name_seq);
      fclose(file_out);
    }
  }
  // Write the unrooted tree
  if(umnode){
    if(NJ){sprintf(nametree, "%s.nj.tree", name_out);}
    else{sprintf(nametree, "%s.al.tree", name_out);}
    printf("Writing unrooted tree in Newick format in %s\n", nametree);
    file_out=fopen(nametree, "w");
    fprintf(file_out, "%s\n", umnode);
    fclose(file_out);
    free(umnode);
  }else{
    printf("WARNING, unrooted tree was not found\n");
  }

  if(umnode_seq){
    char nametree[200];
    if(NJ){sprintf(nametree, "%s.nj.seqtree", name_out);}
    else{sprintf(nametree, "%s.al.seqtree", name_out);}
    printf("Writing unrooted seq tree in Newick format in %s\n", nametree);
    file_out=fopen(nametree, "w");
    fprintf(file_out, "%s\n", umnode_seq);
    fclose(file_out);
    free(umnode_seq);
  }

  for(a=0; a<n_nodes; a++){
    free(div[a]); if(NJ)free(Q[a]);
  }
  if(div_seq){
    for(a=0; a<n_nodes; a++){free(div_seq[a]);}
    //free(div_seq);
  }
  for(a=0; a<n_nodes; a++){
    node=nodes+a;
    free(node->sons);
    free(node->name); 
    if(node->name_seq){free(node->name_seq);}
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

char *Copy_name_3(struct node *node1, struct node *node2, struct node *node3)
{
  char *name=malloc(sizeof(char)*
		    (strlen(node1->name)+
		     strlen(node2->name)+
		     strlen(node3->name)+30));
  sprintf(name, "(%s:%.3g,%s:%.3g,%s:%.3g)",
	  node1->name, node1->branch,
	  node2->name, node2->branch,
	  node3->name, node3->branch);
  return(name);
}

void Copy_name_seq(struct node *mnode, struct node *node1, struct node *node2)
{
  mnode->name_seq=
    malloc((strlen(node1->name_seq)+strlen(node2->name_seq)+60)*sizeof(char));
  sprintf(mnode->name_seq, "(%s:%.3g,%s:%.3g)",
	  node1->name_seq, node1->branch_seq,
	  node2->name_seq, node2->branch_seq);
}

char *Copy_name_3_seq(struct node *node1,struct node *node2,struct node *node3)
{
  char *name=malloc(sizeof(char)*
		    (strlen(node1->name_seq)+
		     strlen(node2->name_seq)+
		     strlen(node3->name_seq)+30));
  sprintf(name, "(%s:%.3g,%s:%.3g,%s:%.3g)",
	  node1->name_seq, node1->branch_seq,
	  node2->name_seq, node2->branch_seq,
	  node3->name_seq, node3->branch);
  return(name);
}

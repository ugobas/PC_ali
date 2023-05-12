/* Listing all maximal cliques
By a result of Moon & Moser (1965), every n-vertex graph has at most 3n/3 maximal cliques. They can be listed by the Bron-Kerbosch algorithm, a recursive backtracking procedure of Bron & Kerbosch (1973). The main recursive subroutine of this procedure has three arguments: a partially constructed (non-maximal) clique, a set of candidate vertices that could be added to the clique, and another set of vertices that should not be added (because doing so would lead to a clique that has already been found). The algorithm tries adding the candidate vertices one by one to the partial clique, making a recursive call for each one. After trying each of these vertices, it moves it to the set of vertices that should not be added again.
However, when the number of cliques is significantly smaller than its worst case, other algorithms might be preferable. As Tsukiyama et al. (1977) showed, it is also possible to list all maximal cliques in a graph in an amount of time that is polynomial per generated clique. An algorithm such as theirs in which the running time depends on the output size is known as an output-sensitive algorithm. Their algorithm is based on the following two observations, relating the maximal cliques of the given graph G to the maximal cliques of a graph G \ v formed by removing an arbitrary vertex v from G:

    For every maximal clique K of G \ v, either K continues to form a maximal clique in G, or K U {v} forms a maximal clique in G. Therefore, G has at least as many maximal cliques as G \ v does.
    Each maximal clique in G that does not contain v is a maximal clique in G \ v, and each maximal clique in G that does contain v can be formed from a maximal clique K in G \ v by adding v and removing the non-neighbors of v from K.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "clique_msa.h"
#include "allocate.h"
#include "rank.h"
#include <limits.h>

float MAX_TEST=2.5;  // The smaller, the faster but less accurate
float MIN_MAX=0.35; // The higher, the faster but less accurate
int joinmax=8;      // Max. attempts to join

int DBG=0;
int Ltot=0;
char **Link=NULL, **Link_indirect=NULL;
int *Prot_site, *Res_site;
int **Site_prot;
int *num_link;
int nprot;
int nmax;
int Cl_max;
int *ini_prot;
int *len_prot;
int print_warning=0;
int *not_stored;
int *max_clique;
int *sum_clique;

struct clique{
  int n;
  int *site;
  int num;
  struct clique *next;
  struct clique *prev;
};

struct rank *rank_clique=NULL, *First_rank=NULL, *Last_rank=NULL;
int **ranked_clique, *num_ranked, max_ranked;

extern void Change_ext(char *, char *, char *);
int Set_sites(int **Prot_site, int **Res_site,
		       int ***Site_prot,
		       int *Lprot, int Nprot);
char **Make_links(int Ltot, int ***Ali_pair, int *Lprot, int Nprot);
char **Make_links_indirect(int Ltot, int ***Ali_pair, int *Lprot, int Nprot);
void Set_link(char **Link, int iprot, int i, int jprot, int j);
int **Make_contacts(int Ltot, char **Link);
void Rank_cc(int *rank_cc, int Ltot, int **Cont_list, char **Link, int Nprot);
void Maximal_cliques(int *n_clique, struct clique *store_clique,
		     int site, int **C_list, char **Link, int Nprot);
int Included(int *sites, int n, struct clique *store_clique);
void Store(int **tmp_store, int *n_store, int n, int cl, int cl_max);
int Find_tmp(int *min_n, int **tmp_store, int *n_store, int n, int cl_max);
struct rank *Set_ranks(int Tot_cliques);

int Order_clique(struct clique *clique,
		 struct clique **First_clique, int Nprot);
int Find_order(struct clique *clique, struct clique *First_clique,
	       struct clique **prev, int *njoin, int *rank_join,
	       int *score, struct clique **cl_join, int joinmax);
int Test_join(struct clique *clique, struct clique *cl);
int Not_assigned(struct clique **First_clique,
		 struct clique *store_clique,
		 int Nprot, int Ltot);
int Test_order(int *ntot, struct clique *First_clique, int *Lprot, int print);
void Write_msa(int **msa, struct clique *First_clique,
	       int Nprot, int *Lprot, int Ltot, int L_msa);
void Print_clique(struct clique *clique);
// Not used:
void Single_sites(int n_clique, struct clique *store_clique,
		  int Ltot, int Nprot);


void Clique_MSA(int **msa, int *L_msa,      // output
		int N_ali_max,
		int ***Ali_pair, // Ali_pair[i][j][site_i]=site_j j<i
		int *Lprot,      // Length of protein i
		int Nprot,       // Number of proteins
		char **Seq,      // All sequences
		char **name_seq, // Names of sequences
		char *name_in,
		int L_msa_ini)   // Length of the input MSA
{
  Ltot=Set_sites(&Prot_site, &Res_site, &Site_prot, Lprot, Nprot);
  nprot=Nprot; len_prot=Lprot; // global variables

  // Make_links
  Link=Make_links(Ltot, Ali_pair, Lprot, Nprot);
  Link_indirect=Make_links_indirect(Ltot, Ali_pair, Lprot, Nprot);
  int **Cont_list=Make_contacts(Ltot, Link);
  //int **Cont_list=Make_contacts(Ltot, Link_indirect);

  // Clustering coefficient
  int rank_cc[Ltot];
  Rank_cc(rank_cc, Ltot, Cont_list, Link, Nprot);

  // max number stored cliques per num. nodes:
  max_ranked=MAX_TEST*Nprot*L_msa_ini;
  Cl_max=max_ranked*Nprot;
  printf("Max number of cliques per n: %d Total: %d\n", max_ranked, Cl_max);
  struct clique store_clique[Cl_max]; int i;
  for(i=0; i<Cl_max; i++)store_clique[i].site=NULL;
  nmax=Nprot+1;
  ranked_clique=malloc(nmax*sizeof(int *));
  num_ranked=malloc(nmax*sizeof(int));
  not_stored=malloc(nmax*sizeof(int));
  for(i=0; i<nmax; i++){
    ranked_clique[i]=malloc(max_ranked*sizeof(int));
    num_ranked[i]=0;
    not_stored[i]=0;
  }
  max_clique=malloc(Ltot*sizeof(int));
  sum_clique=malloc(Ltot*sizeof(int));
  for(i=0; i<Ltot; i++){
    max_clique[i]=0; sum_clique[i]=0;
  }

  printf("Start making cliques for %d sites: ", Ltot); fflush(NULL);
  int n_clique=0;
  for(int i=0; i<Ltot; i++){ 
    int site=rank_cc[i]; // Higher corr coeff first
    if(max_clique[site]>=(num_link[site]*0.5)  ||
       sum_clique[site]>=(num_link[site]*0.75))continue;
    Maximal_cliques(&n_clique, store_clique, site, Cont_list, Link, Nprot);
    //Maximal_cliques(&n_clique,store_clique,site, Cont_list, Link_indirect,
    //		    Nprot);

    printf("%d ", site); fflush(NULL);
  }
  int Tot_cliques=0;
  printf("\n%d cliques ranked: ", n_clique);
  for(i=Nprot; i>1; i--){    
    printf("%d (%d) ", num_ranked[i], i); Tot_cliques+=num_ranked[i];
  }
  printf("  Total=%d\n", Tot_cliques);
  int print_not=1;
  for(i=Nprot; i>1; i--){
    if(not_stored[i]==0)continue;
    if(print_not){printf("not stored cliques: "); print_not=0;}
    printf("%d (%d) ", not_stored[i], i);
  }
  if(print_not==0)printf("\n");

  //Single_sites(n_clique, store_clique, Ltot, Nprot);
  //rank_clique=Set_ranks(Tot_cliques);

  struct clique *First_clique=NULL, *cl;
  int joined=0, ntot=0;
  *L_msa=0;
  //struct rank *r=First_rank;
  //while(r!=NULL){
  //int ncl=r->index;
  for(int n=nprot; n>=1; n--){
    for(i=0; i<num_ranked[n]; i++){
      int ncl=ranked_clique[n][i];
      cl=store_clique+ncl;
      cl->num=Order_clique(cl, &First_clique, Nprot);
      if(cl->num>0){(*L_msa)++; ntot+=cl->n;}
      else if(cl->num==0){joined++; ntot+=cl->n;}
      //r=r->next;
    }
  }
  printf("%d res. %d columns ordered, %d joined initial: %d\n",
	   ntot, *L_msa, joined, L_msa_ini);

  int error=Test_order(&ntot, First_clique, Lprot, 0);
  printf("%d cliques with %d residues set as %d MSA columns, Ltot=%d"
	 " Initial columns= %d\n", *L_msa+joined,ntot,*L_msa,Ltot,L_msa_ini);
  printf("%d errors in Test_order before Not_assigned\n", error);
  if(error){printf("Exiting\n"); exit(8);}

  // Identify sites that have not been assigned to any clique
  (*L_msa)+=Not_assigned(&First_clique,store_clique,Nprot,Ltot);

  error=Test_order(&ntot, First_clique, Lprot, 0);
  printf("%d cliques with %d residues -> %d MSA columns initial: %d Ltot= %d\n",
	 *L_msa+joined, ntot, *L_msa, L_msa_ini, Ltot);
  printf("%d errors in Test_order after Not_assigned\n", error);
  if(error){printf("Exiting\n"); exit(8);}


  if(*L_msa>N_ali_max){
    printf("ERROR, too many columns (%d) in clique_msa. Increase N_ali_max\n",
	   *L_msa); exit(8);
  }
  Write_msa(msa, First_clique, Nprot, Lprot, Ltot, *L_msa);

  char name_msa[100]; Change_ext(name_msa, name_in, ".msa");
  /*printf("Printing clique-based multiple alignment in file %s\n", name_msa);
    Print_MSA(msa, name_msa, Nprot, *L_msa, Lprot, Seq, name_seq); */
  
  // clean
  printf("Cleaning memory in clique_msa\n");
  //printf("Cleaning memory: Sites\n");
  free(Prot_site);
  free(Res_site);
  for(i=0; i<Nprot; i++)free(Site_prot[i]);
  free(Site_prot);
  free(num_link);
  free(max_clique);
  free(sum_clique);

  //printf("Cleaning memory: Links\n");
  int site;
  for(site=0; site<Ltot; site++)free(Link[site]);
  free(Link);
  if(Link_indirect){
    for(site=0; site<Ltot; site++)free(Link_indirect[site]);
    free(Link_indirect);
  }
  //printf("Cleaning memory: Contacts\n");
  for(site=0; site<Ltot; site++)free(Cont_list[site]);
  free(Cont_list);
  //printf("Cleaning memory: cliques\n");
  for(i=0; i<n_clique; i++){ //Cl_max
    free(store_clique[i].site);
  }
  //printf("Cleaning memory: ranks\n");
  for(i=0; i<nmax; i++)free(ranked_clique[i]);
  //printf("Cleaning memory: msa\n");
  //for(i=0; i<Nprot; i++)free(msa[i]);
  //free(msa);

  return;
}

int Set_sites(int **Prot_site, int **Res_site,
		       int ***Site_prot,
		       int *Lprot, int Nprot)
{
  ini_prot=malloc(Nprot*sizeof(int));
  int Ltot=0; int i, j;
  for(i=0; i<Nprot; i++){ini_prot[i]=Ltot; Ltot+=Lprot[i];}

  printf("Total number of residues: %d\n", Ltot);
  if(Ltot>INT_MAX){
    printf("ERROR, > Maximum allowed = %d\n", INT_MAX);
    exit(8);
  }

  *Prot_site=malloc(Ltot*sizeof(int));
  *Res_site= malloc(Ltot*sizeof(int));
  *Site_prot=malloc(Nprot*sizeof(int *));
  int site=0;
  for(i=0; i<Nprot; i++){
    (*Site_prot)[i]=malloc(Lprot[i]*sizeof(int));
    for(j=0; j<Lprot[i]; j++){
      (*Prot_site)[site]=i;
      (*Res_site)[site]=j;
      (*Site_prot)[i][j]=site;
      site++;
    }
  }
  return(Ltot);
}

char **Make_links(int Ltot, int ***Ali_pair, int *Lprot, int Nprot)
{
  char **Link=malloc(Ltot*sizeof(char *));
  int isite, jsite;
  for(isite=0; isite<Ltot; isite++){
    Link[isite]=malloc(Ltot*sizeof(char));
    for(jsite=0; jsite<Ltot; jsite++){
      Link[isite][jsite]='\0';
    }
  }

  for(int iprot=0; iprot<Nprot; iprot++){
    int L=Lprot[iprot];
    for(int jprot=0; jprot<iprot; jprot++){
	int *ali=Ali_pair[iprot][jprot];
	for(int i=0; i<L; i++){
	  int j=ali[i]; if(j<0)continue;
	  Set_link(Link, iprot, i, jprot, j);
	}
    }
  }
  return(Link);
}

void Set_link(char **Link, int iprot, int i, int jprot, int j)
{
  int isite=Site_prot[iprot][i];
  int jsite=Site_prot[jprot][j];
  Link[isite][jsite]='1';
  Link[jsite][isite]='1';
}

char **Make_links_indirect(int Ltot, int ***Ali_pair, int *Lprot, int Nprot)
{
  char **Link=malloc(Ltot*sizeof(char *));
  int isite, jsite;
  for(isite=0; isite<Ltot; isite++){
    Link[isite]=malloc(Ltot*sizeof(char));
    for(jsite=0; jsite<Ltot; jsite++){
      Link[isite][jsite]='\0';
    }
  }

  for(int iprot=0; iprot<Nprot; iprot++){
    int Li=Lprot[iprot], i, j, k;
    for(int jprot=0; jprot<iprot; jprot++){
      int *ali=Ali_pair[iprot][jprot];
      for(int i=0; i<Li; i++){
	int j=ali[i]; if(j<0)continue;
	Set_link(Link, iprot, i, jprot, j);
      }
      for(int kprot=0; kprot<Nprot; kprot++){
	if(kprot==iprot || kprot==jprot)continue;
	if(kprot>iprot){   // jprot<iprot<kprot
	  int *ali_ki=Ali_pair[kprot][iprot];
	  int *ali_kj=Ali_pair[kprot][jprot];
	  for(k=0; k<Lprot[kprot]; k++){
	    i=ali_ki[k]; if(i<0)continue;
	    j=ali_kj[k]; if(j<0)continue;
	    Set_link(Link, iprot, i, jprot, j);
	  }
	}else if(kprot<jprot){ // kprot<jprot<iprot
	    int *ali_ik=Ali_pair[iprot][kprot];
	    int *ali_jk=Ali_pair[jprot][kprot];
	    int Lk=Lprot[kprot], aligned[Lk];
	    for(k=0; k<Lk; k++)aligned[k]=-1;
	    for(i=0; i<Li; i++){
	      k=ali_ik[i]; if(k>=0)aligned[k]=i;
	    }
	    for(j=0; j<Lprot[jprot]; j++){
	      k=ali_jk[j]; if(k<0)continue;
	      i=aligned[k]; if(i<0)continue;
	      Set_link(Link, iprot, i, jprot, j);
	  }
	}else{ // jprot<kprot<iprot
	    int *ali_ik=Ali_pair[iprot][kprot];
	    int *ali_kj=Ali_pair[kprot][jprot];
	    for(i=0; i<Li; i++){
	      k=ali_ik[i]; if(k<0)continue;
	      j=ali_kj[k]; if(j<0)continue;
	      Set_link(Link, iprot, i, jprot, j);
	    }
	}
      }
    }
  }

  return(Link);
}

int **Make_contacts(int Ltot, char **Link)
{
  int tot_link=0;
  int **Cont_list=malloc(Ltot*sizeof(int *)), site;
  num_link=malloc(Ltot*sizeof(int));
  for(site=0; site<Ltot; site++)num_link[site]=0;
  for(site=0; site<Ltot; site++){
    int nlink=0; char *Ls=Link[site]; int s2;
    for(s2=site+1; s2<Ltot; s2++)if(Ls[s2])nlink++;
    Cont_list[site]=malloc((nlink+1)*sizeof(int));
    int *C_list=Cont_list[site];
    nlink=0;
    for(s2=site+1; s2<Ltot; s2++){
      if(Ls[s2]){C_list[nlink]=s2; num_link[s2]++; nlink++;}
    }
    C_list[nlink]=-1;
    tot_link+=nlink;
    num_link[site]+=nlink;
  }
  printf("Mean number of aligned seq per site: %.0f (max= %d)\n",
	 2.0*(float)tot_link/Ltot, nprot-1);
  if(0){
    for(site=0; site<20; site++){
      printf("%d aligned to site %d: ",num_link[site], site);
      int *s2=Cont_list[site];
      while(*s2>=0){printf("%d ",*s2); s2++;}
      printf("\n");
    }
    exit(8);
  }

  return(Cont_list);
}
    
void Maximal_cliques(int *n_clique, struct clique *store_clique,
		     int site, int **C_list, char **Link, int Nprot)
{
  int *s2=C_list[site];
  if(*s2<0)return; // No links are present
  int cl_max=2*max_ranked, warning=0;
  int nmax=Nprot+1;

  struct clique tmp_clique[cl_max];
  int *tmp_store[nmax], n_store[nmax], i, j;
  for(i=0; i<nmax; i++){
    tmp_store[i]=malloc(cl_max*sizeof(int));
    n_store[i]=0;
  }
  int min_n=1, *nodes[cl_max], n_nodes[cl_max];
  for(i=0; i<cl_max; i++)nodes[i]=NULL;

  // Store site in clique cl=0
  int cl=0;
  tmp_clique[cl].n=1;
  tmp_clique[cl].site=malloc(Nprot*sizeof(int));
  tmp_clique[cl].site[0]=site;
  Store(tmp_store, n_store, 1, cl, cl_max);
  cl++;
  
  while(*s2>=0){

    int cl_new=cl;
    if(max_clique[*s2]>=num_link[*s2]/2  ||
       sum_clique[*s2]>=num_link[*s2]*0.75)goto next_link;

    // Loop on all previously stored cliques that do not include *s2
    int added=0;
    for(i=0; i<cl; i++){ 
      struct clique *cli=tmp_clique+i;
      int m=0, curr_nodes[Nprot];
      for(j=0; j<cli->n; j++){
	if(Link[cli->site[j]][*s2]){ // Link_indirect
	  curr_nodes[m]=cli->site[j]; m++;
	}
      }
      curr_nodes[m]=*s2;
      /*if(site>=len_prot[0] && Included(curr_nodes, m+1, store_clique)){
	continue; // larger clique has been already stored
	}*/
      if(m==cli->n){
	cli->site[m]=*s2; cli->n++; // all linked, add *s2 to clique cli
      }else if(m <= min_n){
	continue;  
      }else{          // only m<n linked
	// Check if already used
	int exists=0;
	for(int k=0; k<added; k++){
	  if(n_nodes[k]<m)continue;
	  int *nk=nodes[k];
	  for(j=0; j<m; j++){
	    if(nk[j]!=curr_nodes[j])break;
	  }
	  if(j==m){exists=1; break;}
	}
	if(exists)continue;
	
	// Make another clique
	struct clique *cl2; int n=m+1;
	if(cl_new==cl_max){
	  // Write clique at the place k of the minimum previous clique
	  // that is < n (if there is not, k=-1)
	  if(print_warning==0){
	    printf("WARNING, too many candidate cliques > %d for site %d. ",
		   cl_max, site); print_warning=1;
	  }
	  warning ++;
	  int k=Find_tmp(&min_n, tmp_store, n_store, n, cl_max);
	  if(k<0){continue;}
	  cl2=tmp_clique+k;
	  Store(tmp_store, n_store, n, k, cl_max);
	}else{
	  cl2=tmp_clique+cl_new;
	  Store(tmp_store, n_store, n, cl_new, cl_max);
	  cl2->site=malloc(Nprot*sizeof(int));
	  cl_new++;
	}
	for(j=0; j<m; j++)cl2->site[j]=curr_nodes[j];
	cl2->site[m]=*s2;
	cl2->n=n;
	// store
	if(added<cl_max){
	  if(nodes[added]==NULL)nodes[added]=malloc(Nprot*sizeof(int));
	  for(j=0; j<m; j++){nodes[added][j]=curr_nodes[j];}
	  n_nodes[added]=m;
	  added++;
	}
      } // end make clique
    } // end loop on previous cliques
  next_link:
    s2++;
    cl=cl_new;
  }
  if(warning){
    printf("%d cliques not stored %d stored (%.4f) for site %d ",
	   warning, cl, cl/(float)(warning+cl), site);
  }
  // Clean memory
  for(i=0; i<cl_max; i++){if(nodes[i])free(nodes[i]);}

  // Omit cliques with fewer than half the maximum nodes 
  int max_n=0;
  for(i=0; i<cl; i++)if(tmp_clique[i].n>max_n)max_n=tmp_clique[i].n;
  min_n=max_n*MIN_MAX;
  
  for(i=0; i<cl; i++){
    struct clique *cli=tmp_clique+i;
    if(cli->n <= min_n)continue;
    if(nprot<40 && site>=len_prot[0] && //not checked if nprot>40
       Included(cli->site, cli->n, store_clique)){
      continue; 
    }

    // k: where to store
    // k=Rank(cli->n,rank_clique,n_clique,Cl_max,First_rank,Last_rank);
    //if(k<0)continue;
    if(num_ranked[cli->n]>=max_ranked){
      not_stored[cli->n]++; continue;
    }
    int k=*n_clique; (*n_clique)++;
    if(k>Cl_max){
      printf("ERROR in Make_cliques, too many cliques (> %d) n= %d\n",
	     Cl_max, cli->n); exit(8);
    }
    ranked_clique[cli->n][num_ranked[cli->n]]=k;
    num_ranked[cli->n]++;
    
    // store clique at k
    struct clique *store=store_clique+k;
    if(store->site==NULL)store->site=malloc(nprot*sizeof(int));
    for(j=0; j<Nprot; j++)store->site[j]=-1;
    for(j=0; j<cli->n; j++){
      int site=cli->site[j];
      int iprot=Prot_site[site];
      store->site[iprot]=site;
      if(cli->n>max_clique[site])max_clique[site]=cli->n;
      sum_clique[site]+=(cli->n-1);
    }
    store->n=cli->n;
    store->num=1;
  }
  //printf(" tot_cliques=%d\n", *n_clique);
  
  for(i=0; i<cl; i++)free(tmp_clique[i].site);
  for(i=0; i<nmax; i++)free(tmp_store[i]);
}

int Included(int *site, int n, struct clique *store_clique)
{
  // Check if click i (still in tmp format) is included in one of the
  // previous clicks with more nodes
  int iprot[n], j;
  for(j=0; j<n; j++){iprot[j]=Prot_site[site[j]];}
  for(int nn=n+1; nn<=nprot; nn++){
    for(int k=0; k<num_ranked[nn]; k++){
      struct clique *cl=store_clique+ranked_clique[nn][k];
      for(j=0; j<n; j++){
	if(site[j]!=cl->site[iprot[j]])break;
      }
      if(j==n)return(1); // Clique cli is included in cl
    }
  }
  return(0);  // Clique cli is not included in any cl
}

void Store(int **tmp_store, int *n_store, int n, int cl, int cl_max)
{
  if(n<=0 || n>nprot){
    printf("ERROR, wrong n= %d in clique %d over %d\n", n, cl, cl_max);
    exit(8);
  }
  tmp_store[n][n_store[n]]=cl; n_store[n]++;
  if(n_store[n] > cl_max){
    printf("ERROR, too many cliques with %d nodes: %d > %d (last: %d)\n",
	   n, n_store[n], cl_max, cl); exit(8);
  }
}

int Find_tmp(int *min_n, int **tmp_store, int *n_store, int n, int cl_max)
{
  int nn; // find lowest nn with stored cliques
  for(nn=2; nn<n; nn++)if(n_store[nn])break;
  *min_n=nn;
  if(n<nn)return(-1);
  int k=n_store[nn]-1; if(k<0)return(-1);
  tmp_store[nn][k]=-1; n_store[nn]--;  // eliminate from nn
  //tmp_store[n][n_store[n]]=k; n_store[n]++; // add to n
  return(k);
}

void Single_sites(int n_clique, struct clique *store_clique,
		  int Ltot, int Nprot)
{
  int cl_tot=n_clique+Ltot;
  if(cl_tot > Cl_max){
    printf("ERROR in Single sites, too many cliques\n"); exit(8);
  }
  
  ranked_clique[1]=malloc(Ltot*sizeof(int));
  num_ranked[1]=Ltot;
  for(int site=0; site<Ltot; site++){
    int i=site+n_clique;
    ranked_clique[1][site]=i;
    struct clique *cli=store_clique+i;
    cli->n=1;
    cli->site=malloc(Nprot*sizeof(int));
    for(int j=0; j<Nprot; j++)cli->site[j]=-1;
    cli->site[Prot_site[site]]=site;
  }
}

struct rank *Set_ranks(int Tot_cliques){
  struct rank *ranks=malloc(Tot_cliques*sizeof(struct rank)), *r=ranks;
  First_rank=r;
  for(int n=nprot; n>=1; n--){
    for(int i=0; i<num_ranked[n]; i++){
      r->score=n;
      r->index=ranked_clique[n][i];
      r->next=r+1;
      r++;
    }
  }
  r--; Last_rank=r; r->next=NULL;
  return(ranks);
}

int Order_clique(struct clique *clique, struct clique **First_clique,
		 int Nprot)
{
  // 1: ordered 0: joined -1: unassigned -2: not allowed

  if(*First_clique==NULL){
    *First_clique=clique;
    clique->next=NULL;
    clique->prev=NULL;
    return(1);
  }
  
  int score[joinmax], rank_join[joinmax], njoin=0, j;
  struct clique *prev=NULL, *cl_join[joinmax];
  if(DBG)printf("Finding order\n"); //DBG
  // njoin<=joinmax is the
  int order=Find_order(clique, *First_clique, &prev, 
		       &njoin, rank_join, score, cl_join, joinmax);
  if(DBG)printf("Order= %d\n", order); //DBG
  if(order<0)return(order);
  
  if(njoin <= 0)goto link;

  // Join the cliques clique and cl_join
  if(DBG)printf("Trying to join with %d candidates\n", njoin); //DBG
  struct clique tmp, *cl; //, *first=*First_clique; 
  tmp.site=malloc(Nprot*sizeof(int));

  int kmax=-1;
  for(int i=0; i<njoin; i++){
    int k=rank_join[i];
    if(k<0 || k>= Cl_max)continue;
    if(score[k]<0 || score[k]>=Nprot)continue;
    cl=cl_join[k];
    if(cl==NULL || cl->n <= 0 || cl->n > Nprot)continue;
    tmp.n=cl->n+clique->n;
    for(j=0; j<Nprot; j++){
      if(clique->site[j]>=0){tmp.site[j]=clique->site[j];}
      else{tmp.site[j]=cl->site[j];}
    }
    order=Find_order(&tmp,*First_clique,NULL,NULL,NULL,NULL,NULL,0);
    if(order>0){kmax=k; break;} // accept

    if(DBG)printf("Testing joined clique\n");  //DBG
    /*tmp.prev=cl->prev;
    tmp.next=cl->next;
    if(tmp.prev){tmp.prev->next=&tmp;}
    else{*First_clique=&tmp;}*/
    int error=0; //Test_order(NULL, *First_clique, len_prot, 0);
    /*if(error){printf("ERROR in joined clique\n");}
    if(DBG)printf("Done\n");  //DBG
    if(tmp.prev){tmp.prev->next=cl;}
    else{*First_clique=cl;}*/
    if(error==0){kmax=k; break;}
  }
  free(tmp.site);
  if(order==-3 || kmax<0)goto link;

  cl=cl_join[kmax];
  cl->n += clique->n;
  if(cl->n>Nprot){
    printf("ERROR, too many proteins in joined clique: %d > %d\n",
	   cl->n, Nprot); exit(8);
  }
  for(j=0; j<Nprot; j++){
    if(clique->site[j]>=0)cl->site[j]=clique->site[j];
  }
  cl->num++;
  if(DBG)printf("Joined\n");  //DBG
  return(0);

 link:
  if(DBG)printf("Not joined, linking\n");  //DBG
  if(prev==NULL){
    clique->next=*First_clique;
    *First_clique=clique;
  }else{
    clique->next=prev->next;
    prev->next=clique;
  }
  clique->prev=prev;
  if(clique->next)clique->next->prev=clique;
  if(DBG)printf("Linked\n");  //DBG

  return(1);

}

int Not_assigned(struct clique **First_clique,
		 struct clique *store_clique,
		 int Nprot, int Ltot)
{
  // Find unassigned cliques for storing unassigned sites
  int error=0, i;

  printf("Ordering not assigned sites\n");
  int site;
  int assigned[Ltot]; for(site=0; site<Ltot; site++)assigned[site]=0;
  struct clique *clique=*First_clique;
  while(clique != NULL){
    for(i=0; i<Nprot; i++){
      site=clique->site[i];
      if(site >= Ltot){
	printf("ERROR, site %d > %d]\n", site+1, Ltot); error=1;
      }
      if(site >=0)assigned[site]=1;
    }
    clique=clique->next;
  }

  if(DBG)printf("Counting not assigned sites\n");

  // Count not assigned, rank them based on number of links
  int not_assigned=0;
  int *ranked[nmax], rmax=Ltot;
  for(i=0; i<nmax; i++){
    num_ranked[i]=0;
    ranked[i]=malloc(rmax*sizeof(int));
  }
  for(site=Ltot-1; site>=0; site--){
    if(assigned[site])continue;
    not_assigned++;
    int i=num_link[site]; if(i>nprot)i=nprot;
    ranked[i][num_ranked[i]]=site;
    num_ranked[i]++;
    if(num_ranked[i]>=rmax){
      printf("ERROR, > %d not assigned sites with %d links\n",
	     rmax, i); exit(8);
    }
  }
  if(DBG)printf("%d not assigned sites\n", not_assigned);

  // Prepare where to store
  int n_empty=0, empty[Cl_max];
  for(i=0; i<Cl_max; i++){
    if(store_clique[i].num<=0){
      empty[n_empty]=i; n_empty++;
    }
  }
  if(DBG)printf("%d empty boxes for storing them\n", n_empty);

  int added=0, joined=0;
  for(i=nprot; i>=0; i--){ // From more to fewer links 
    //for(i=0; i<=nprot; i++){ // From fewer to more links
    for(int ii=0; ii<num_ranked[i]; ii++){
      site=ranked[i][ii];
      if(site<0 || site>=Ltot){
	printf("ERROR, wrong site %d %d links ii=%d over %d\n",
	       site, i, ii, num_ranked[i]); exit(8);
      }
      if(assigned[site]){
	printf("WARNING, site %d is already assigned\n", site);
	continue;
      }
      int k=empty[added];
      if(added >= n_empty || k<0 || k>= n_empty){
	printf("ERROR, not assigned site %d %d %d cannot be stored at %d "
	       "not in [0,%d] n_empty= %d\n",
	       not_assigned, site, added, k, Cl_max, n_empty);
	error=1;
      }
      //printf("k=%d\n", k); fflush(NULL);
      struct clique *cl=store_clique+k;
      cl->n=1;
      if(cl->site==NULL)cl->site=malloc(nmax*sizeof(int));
      for(int j=0; j<Nprot; j++)cl->site[j]=-1;
      cl->site[Prot_site[site]]=site;
      cl->num=Order_clique(cl, First_clique, Nprot);
      if(cl->num<0){
	printf("ERROR, not allowed site %d prot=%d res=%d\n",
	       site, Prot_site[site], Res_site[site]);
	error=1;
      }else if(cl->num>0){
	added++;
	if(added > n_empty){
	  printf("ERROR, too few empty cliques (%d) for storing "
		 "not assigned sites\n", n_empty); error=1;
	}
      }else{
	joined++;
      }
    }
  }
  for(i=0; i<nmax; i++)free(ranked[i]);

  printf("%d unassigned sites, of which %d inserted columns %d joined\n",
	 not_assigned, added, joined);
  if(error){printf("Exiting because of errors\n"); exit(8);}
  return(added);

}

int Test_order(int *ntot, struct clique *First_clique, int *Lprot, int print)
{
  int error=0, i, k=0; if(ntot)*ntot=0;
  struct clique *cl=First_clique; 
  int prev_res[nprot]; for(i=0; i<nprot; i++)prev_res[i]=-1;
  while(cl!=NULL){
    if(print){Print_clique(cl);}
    for(i=0; i<nprot; i++){
      int r=cl->site[i]; if(r<0)continue;
      if(r<prev_res[i]){
	printf("ERROR cl %d prot %d L=%d, res %d located after %d n=%d\n",
	       k, i, Lprot[i], r-ini_prot[i], prev_res[i]-ini_prot[i], cl->n);
	error++;
      }
      prev_res[i]=r;
    }
    if(ntot){(*ntot)+=cl->n;} cl=cl->next; k++;
  }
  return(error);
}

void Write_msa(int **msa, struct clique *First_clique,
	       int Nprot, int *Lprot, int Ltot, int L_msa)
{
  int error=0, i;
  for(i=0; i<Nprot; i++){
    for(int j=0; j<L_msa; j++)msa[i][j]=-1;
  }
  struct clique *cl=First_clique; int col=0;
  while(cl!=NULL){
    for(i=0; i<Nprot; i++){
      int site=cl->site[i];
      if(site>=0){
	if(site>= Ltot){
	  printf("ERROR in Write_msa, site %d not in [0,%d]\n",
		 site, Ltot-1); error=1;
	}
	if(Prot_site[site]!= i){
	  printf("ERROR in Write_msa, protein %d != %d\n",
		 Prot_site[site], i); error=1;
	}
	int res=Res_site[site];
	if(res<0 || res>= Lprot[i]){
	  printf("ERROR in Write_msa, residue %d not in [0,%d]\n",
		 res, Lprot[i]-1); error=1;
	}
	msa[i][col]=res;
      }else{
	msa[i][col]=-1;
      }
    }
    cl=cl->next; col++;
  }
  if(col != L_msa){
    printf("ERROR in Write_msa, found %d columns instead of %d\n",
	   col, L_msa); error=1;
  }
  if(error){printf("Exiting\n"); exit(8);}
}


void Print_MSA(int **msa_new, char *name_msa, int Nprot, int L_msa,
	       int *Lprot, char **Seq, char **name_seq)
{
  printf("Printing MSA in file %s\n",name_msa);
  FILE *file_out=fopen(name_msa, "w");
  for(int iprot=0; iprot<Nprot; iprot++){
    fprintf(file_out, ">%s\n", name_seq[iprot]);
    int *ali=msa_new[iprot];
    for(int i=0; i<L_msa; i++){
      if(ali[i]>=0){
	fprintf(file_out, "%c", Seq[iprot][ali[i]]);
      }else{
	fprintf(file_out, "-");
      }
    }
    fprintf(file_out, "\n");
  }
  fclose(file_out);
}
void Print_clique(struct clique *clique){
  printf("%2d  ", clique->n);
  for(int j=0; j<nprot; j++){
    int res=clique->site[j];
    if(res>=0){printf(" %3d", res-ini_prot[j]);}
    else{printf("   -");}
  }
  printf(" c=%d", clique->num);
  printf("\n");
}

int Find_order(struct clique *clique, struct clique *First_clique,
	       struct clique **prev_cl,
	       int *njoin, int *rank_join, int *score,
	       struct clique **cl_join, int joinmax)
{
  // 1: ok, can be joined or ordered -2: repeated -3: wrong order
  if(njoin)*njoin=0;

  int left=0; if(prev_cl){*prev_cl=NULL;} 
  int n_left=0, n_right=0, n_overlap=0, n_tot=0, j;
  int n_l=0, n_r=0, n_ov=0, tot=0;
  int join1=joinmax-1;

  struct clique *cl=First_clique, *prev=NULL;
  while(cl != NULL){
    
    // Compare proteins in clique (target) and cl (ordered)
    //int n_downstream=0, n_upstream=0, n_overlap=0, tot=0;
    n_l=0, n_r=0; n_ov=0; tot=0;
    for(j=0; j<nprot; j++){
      if(clique->site[j]>=0 && cl->site[j]>=0){
	tot++;
	if(clique->site[j]<cl->site[j]){n_l++;}
	else if(clique->site[j]>cl->site[j]){n_r++;}
	else{n_ov++;}
      }
    }
    
    if(n_ov && joinmax)return(-2); // site already assigned, discard
    if(tot==0){
      /* No shared proteins, test if the cliques can be joined
	 If all proteins are different and the previous res of
	 the protein is < res and there are > 1 link, join cliques
      */
      if(left || joinmax==0){goto next_cl;}

      // More relaxed conditions for only one res: indirect links, 1 link
      char **links; int min_m;
      if(clique->n==1){links=Link_indirect; min_m=0;}  // Link_indirect
      else{links=Link_indirect; min_m=1;}
      
      // Try to join clusters
      // Compute score:
      if(DBG)printf("Trying to join\n");  //DBG
      int m=0; // number of links
      for(j=0; j<nprot; j++){
	int s1=clique->site[j]; if(s1<0)continue;
	for(int k=0; k<nprot; k++){
	  if(k==j)continue;
	  int s2=cl->site[k];
	  if(s2>=0 && (links[s1][s2]))m++;
	}
      }
      if(m<min_m)goto next_cl; // not joined
      // n.links m>=min_m: find box (j) and rank (k)
      // clique with rank join1: cl=j score=num. links
      if(*njoin >= joinmax){
	if(m <= score[rank_join[join1]])goto next_cl;
	j=rank_join[join1];
      }else{
	j=*njoin;
      }
      int k=0, l;
      for(k=0; k<(*njoin); k++){
	if(m > score[rank_join[k]])break;
      }
      if(k >= joinmax)goto next_cl;
      // Test order of joined clique
      if(DBG)printf("Testing join\n");  //DBG
      if(Test_join(clique, cl)==0)goto next_cl;
      if(DBG)printf("Tested\n");  //DBG

      // store in box j with rank k
      score[j]=m; cl_join[j]=cl; 
      for(l=*njoin; l>k; l--)rank_join[l]=rank_join[l-1];
      rank_join[k]=j;
      if(*njoin < joinmax){(*njoin)++;}

      // end join
      
    }else if(n_l){
      if(n_l!=tot){return(-3);}// Ambigous situation
      if(left==0){
	if(prev_cl){*prev_cl=prev;} left=1;
	n_tot=tot; n_left=n_l; n_right=n_r; n_overlap=n_ov; 
      }
      if(n_l==clique->n)break;
    }else if(n_r){ 
      if(left || n_l){return(-3);}// Ambigous situation
    }
  next_cl:
    prev=cl; cl=cl->next;
  }
  if(left==0){if(prev_cl){*prev_cl=prev;}}

  if(0){
    if(left==0){
      n_tot=tot; n_left=n_l; n_right=n_r; n_overlap=n_ov;
    }
    printf("%d n_left=%d n_right=%d overlap=%d\n",
	   n_tot, n_left, n_right, n_overlap);
    if(prev_cl){
      if(*prev_cl){printf("prev: %d", (*prev_cl)->n);}
      else{printf("prev: NULL");}
      if(*prev_cl)Print_clique(*prev_cl);
    }
    Print_clique(clique);
  }

  return(1);
}

int Test_join(struct clique *clique, struct clique *cl)
{
  int tested[nprot], j;
  // test to the right 
  struct clique *cl2=cl->next;
  for(j=0; j<nprot; j++)tested[j]=0;
  while(cl2 !=NULL){
    int t=0;
    for(j=0; j<nprot; j++){
      if(clique->site[j]<0 || tested[j])continue;
      if(cl2->site[j]>=0){
	if(cl2->site[j]<=clique->site[j]){return(0);} // order fails
	tested[j]=1;
      }
      t++;
    }
    if(t==0)break; // nothing more to test
    cl2=cl2->next;
  }

  // test to the left 
  cl2=cl->prev;
  for(j=0; j<nprot; j++)tested[j]=0;
  while(cl2 !=NULL){
    int t=0;
    for(j=0; j<nprot; j++){
      if(clique->site[j]<0 || tested[j])continue;
      if(cl2->site[j]>=0){
	if(cl2->site[j]>=clique->site[j]){return(0);} // order fails
	tested[j]=1;
      }
      t++;
    }
    if(t==0)break; // nothing more to test
    cl2=cl2->prev;
  }
  return(1);
}

void Rank_cc(int *rank_cc, int Ltot, int **Cont_list, char **Link, int Nprot)
{
  int site, num_cc[Ltot];
  for(site=0; site<Ltot; site++)num_cc[site]=0;
  for(site=0; site<Ltot; site++){
    int *s1=Cont_list[site];
    while(*s1>=0){
      int *s2=s1+1;
      while(*s2>=0){
	if(Link[*s1][*s2]!='\0'){
	  num_cc[site]++; num_cc[*s1]++; num_cc[*s2]++;
	}
	s2++;
      }
      s1++;
    }
  }
  int ccmax=0;
  for(site=0; site<Ltot; site++)if(num_cc[site]>ccmax)ccmax=num_cc[site];
  int rank[ccmax+1][Ltot], nrank[ccmax+1], ranked=0, n;
  for(n=0; n<=ccmax; n++)nrank[n]=0;
  for(site=0; site<Ltot; site++){
    n=num_cc[site]; rank[n][nrank[n]]=site; nrank[n]++;
  }
  for(n=ccmax; n>=0; n--){
    for(int i=0; i<nrank[n]; i++){rank_cc[ranked]=rank[n][i]; ranked++;}
  }
  printf("%d sites ranked by corr coeff, maximum= %d\n", ranked, ccmax);
  if(ranked!=Ltot){
    printf("ERROR in Rank_cc, %d sites but %d ranked\n", Ltot, ranked);
    exit(8);
  }
}

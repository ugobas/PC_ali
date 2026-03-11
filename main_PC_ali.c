/* 
   Program PC_ali
   Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa (CSIC-UAM)
   ubastolla@cbm.csic.es
   PC_ali performs hybrid multiple structure and sequence alignments based on
   the structure+sequence similarity score PC_sim, prints pairwise similarity
   scores and divergence scores and neighbor-joining phylogenetic tree obtained
   with the hybrid evolutionary divergence based on PC_sim.
   Optionally, it computes violations of the molecular clock for each pair
   of proteins.

   It takes as input either not aligned sequences (option -seq) or an MSA
   (option -ali). PDB file names must be specified as sequence names.
   Performs hybrid multiple structure and sequence alignment and computes
   neighbor joining tree based on the hybrid PC_Div divergence.

   Usage:
   PC_ali -seq <sequences in FASTA format, with names of PDB files>
   -ali <MSA file in FASTA format, with names of PDB files>
   # The pdb code is optionally followed by the chain index
   # Ex: >1opd.pdb A or >1opdA or >1opd_A
   -pbdir <path of pdb files>  (default: current folder)
   -pdbext <extension of pdb files>  (default: .pdb)
   -pdblist <File with PDB codes, chains and domain decomposition>

   OUTPUT: MSA (.msa), NJ tree (.tree), structural similarity scores (.sim)
   and structural divergence (.div) for each protein pair, correlations between
   different types of sequence and structure identity (.id), MSA of secondary
   structure (.ss.msa)

*/

#include "Contact_divergence_aux.h"
#include "D_Cont.h"
#include "protein.h"
#include "cont_list.h"
#include "allocate.h"
//#include "normalization.h"
#include "tm_score.h"
#include "read_structures.h"
#include "tree.h"
#include "CV_statistics.h"
#include "align_ss.h"
#include "PC_ali.h"
//#include "consensus_msa.h"
#include "clique_msa.h"
#include "Print_pairwise.h"
#include "nj_align.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
double nbtops = CLOCKS_PER_SEC;
float seq_thr=0.95;  // join sequences with > seq_thr identity

int DIV_ALIGNED=0; // Use aligned fraction score for computing divergence?
int MAKE_CLIQUE=0; // Make clique as initial MSA?
int OMIT_SAME=0;   // Consider only first conformation for same seq group
int PRINT_CV=0;   // Print clock violations?
int PRINT_PDB=0;  // Print structures?
int PRINT_ID=0;  // Print statistics of conserved residues?
//int PRINT_PAIR=0; // Print pairwise alignments?
//int PRINT_CLIQUE=1; // Print multiple alignments based on cliques of pairwise?
int ALI_SS=0;     // Modify Input alignment considering sec.str.?
int ALI_CO=1;     // Modify Input alignment targetting CO?
int ALI_TM=1;     // Modify Input alignment targetting TM?
int ALI_PC=1;     // Modify Input alignment targetting PC?
int PC_NW=1;      // Align PC with NW instead of closest neighbors 
int ITMAX=3;     // Rounds of optimization of target similarity scores
int NMSA=9;     // How many MSA are constructed
int NMSA_act;
#define NTYPE 16     // > 6+NMSA;
int AV_LINK=1;   // Make tree with average linkage (1) of neighbor joining (0)
int CLIQUE_SS=0; // Correct cliques with secondary structure? Not good...
int INP_MSA=0;   // Sequences input as MSA?

int PRINT_SIM=0;  // Print matrix of similarities?
int PRINT_DIV=0;
int PRINT_AVE=0;  // Print average values as a function of d?
int PRINT_GAP=0; // Examine relationship between Triangle Inequality and gaps?
int PRINT_SIM_ALL=0; // Print similarity for all conformations?
int ALL_PAIRS=0; // write all pairwise alis or only j<i?
int ALL_CONF=1;  // Compute struct sim based on all conformations

int PRINT_FASTA=0; // Print PDB sequencs of selected domains
int EXIT_FASTA=0;  // Exit after printing FASTA file

float S0=0.05;    // for Tajima-Nei divergence
float TM0=0.167; //  TM score of unrelated proteins

//float alpha=0.7; // Exponent for computing clock violations


// Parameters for secondary structure correction
int SS_MULT=1;    // sec.str. targeted multiple or pairwise alignment
int SHIFT_MAX=9;  // Maximum allowed shift when targeting sec.str.
float DMAX=10;    // Very large divergence
//float TM0=0.17;   // TM score of unrelated proteins
int NORM, NORMA;

#define EXT_DIFF ".diff" // Extension for comparison of alignments
#define EXT_DIV ".div"   // Extension for divergence output file
#define EXT_SIM ".sim"   // Extension for similarity output file
#define EXT_CV  ".cv"    // Extension for clock-violation output file

#define IJ_MIN_DEF 3
#define CONT_TYPE_DEF 'c'
#define CONT_THR_DEF 4.5

float fun_high=0.90;
float fun_low=0.70;

int IJ_MIN=IJ_MIN_DEF;        // Only contacts with |i-j|>=IJ_MIN
float CONT_THR=CONT_THR_DEF;
char CONT_TYPE=CONT_TYPE_DEF;  //a=alpha b=beta c=all atoms
char CONT_STRING[80];
char CODENAME[40]="main_PC_ali.c";
int CONT_DEF=1;
char pdblist[80]="";

struct Prot_input{
  char code[80];
  char domname[40];
  char domain[40];
  char chain;
  int chain_num;
  char *seq;
  int len;
  int dir;
  int nfrag;
  int *ini_frag;
  int *end_frag;
};
int CHAIN_NUM=0;

// Alignments
int ALL_TYPES=0; // Record for all ali types or only INP and PCAli-mult?
int ATYPE=0, PTYPE=0, CTYPE=0, PT1;
char *ali_name[NTYPE+1];
int code[NTYPE];
static void Set_type(char **ali_name,int *code,int *ATYPE,char *name,int it);
extern char *SS_code;

// Statistics
#define IBIN 10
#define IBIN1 11
int ini_sum=0;
double sum_SS_SI[2], sum_SS_TM[2], sum_SS_CO[2];
double *sum_SS_SI_TM[2], *sum_SS_SI_CO[2], *sum_SS_TM_CO[2];
double // aligned:
**sum_all_ct[NTYPE], //[2][IBIN1],
  **sum_ali_ct[NTYPE], //[2][IBIN1],
  **sum_aaid_ct[NTYPE], //[2][IBIN1],
  **sum_sup_ct[NTYPE], //[2][IBIN1],
  **sum_aaid_sup_ct[NTYPE], //[2][IBIN1],
  **sum_noid_nosup_ct[NTYPE], //[2][IBIN1],
  **sum_ali_cont_ct[NTYPE], //[2][IBIN1],
  **sum_id_cont_ct[NTYPE], //[2][IBIN1],
  // Neighbors:
  **sum_neigh_ali_ct[NTYPE], //[2][IBIN1],
  **sum_neigh_noali_ct[NTYPE], //[2][IBIN1],
  **sum_neigh_noali_aaid_ct[NTYPE], //[2][IBIN1],
  **sum_sup_noneigh_ct[NTYPE], //[2][IBIN1],
  **sum_shift_sup_noneigh_ct[NTYPE], //[2][IBIN1],
  // contacts
  **sum_ali_cont_sup[NTYPE], //[3][IBIN1],
  **sum_id_cont_sup[NTYPE], //[3][IBIN1],
  **sum_ali_cont_aaid[NTYPE], //[3][IBIN1],
  **sum_id_cont_aaid[NTYPE], //[3][IBIN1],
  **sum_ali_cont_aaid_sup[NTYPE][3], //[2][IBIN1],
  **sum_id_cont_aaid_sup[NTYPE][3];  //[2][IBIN1];

double c_ave; // Mean number of contacts per residue

//double PC_load[4]={0.83,0.80,0.94,0.95}; // ali, SI, TM, CO
double PC_load[4]={0.80,0.85,0.92,0.97};
double PC0_1, PC_norm, PC_norm_div;
int opt_PC[NTYPE], Diff_opt[NTYPE];
int npair_pdb=0, npair_seq=0;
double sum_norm_c;
int NSIM=5; // nali SI TM CO PC
double *Sim_ave[5], **Sim_diff[5], **Sim_diff_2[5];

int *rep_str=NULL;

int *id_aa[NTYPE], *id_sup[NTYPE], *shift[NTYPE],
  *neigh_ali[NTYPE], *neigh_noali[NTYPE], *neigh_noali_aaid[NTYPE],
  *ali_cont_ct[NTYPE], *id_cont_ct[NTYPE], 
  *ali_cont_sup[NTYPE], *id_cont_sup[NTYPE],
  *ali_cont_aaid[NTYPE], *id_cont_aaid[NTYPE];


extern int Cluster_seqs(int **rep_str, // repr. str. (one per cluster)
			int **N_all, // total conf per clus
			int ***rep_conf, int ***rep_index, int **N_conf,
			// repr. conf (N_conf per clus)
			int ****all_conf, int ***N_conf_all,
			// (N_conf per clus)
			int ****ali_clus, int **L_ali, // Mult. alignments
			int ****ali_pairs,
			float ****TM_pairs, // TM scores of cluster i 
			struct protein *prots, int N_pdb, // input
			float seq_thr, float str_thr);
extern void Print_pdb(char *name_in, struct protein **prot_p, int n);
extern float Compute_RMSD(struct protein *pi, struct protein *pj,
			  int si, int sj, int **msa, int L_msa);
extern void All_distances(float **d2, float *x1, int n1, float *x2, int n2);

void Score_alignment(float *nali, float *SI, float *TM, float *CO, float *PC,
		     int i, int j, float **Seq_diff,
		     float **na_all, float **SI_all, float **TM_all,
		     float **CO_all, float **PC_all,
		     float **Rot, float *Shift,
		     float c_ave, float *d02, float **d2,
		     int *ali_all, int *ali_ij, int Comp_TM, float TMs,
		     struct protein *proti, struct protein *protj,
		     float norm_c, float norm_ali, int it, int it_opt);

void Sum_scores(double *nali_sum, double *SI_sum, double *TM_sum,
		double *CO_sum, double *PC_sum,
		float *nali, float *SI, float *TM,
		float *CO, float *PC, int it);
void help(char *pname);
void Get_input(char *file_ali, char *file_list, char *file_fun,
	       char *name, char *PDBDIR, char *PDBEXT, char *OUTG,
	       int *NORM, int *ALI_SS, int *SHIFT_MAX, int *SS_MULT,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV,
	       //int *PRINT_PAIR, int *PRINT_CLIQUE,
	       int argc, char **argv);
int Get_sequences(struct Prot_input **Prot_input, int *Nali, char *file_ali,
		  int INP_MSA); //char *PDBDIR, char *PDBEXT, 
int Get_pdb_list(struct Prot_input **Prot_input, char *input);
char *Name_nodir(char *file);
void Name_noext(char *file_new, char *file);
int Read_domain(int **ini_frag, int **end_frag, char *domain);
int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N);
int Remove_gap_cols(int **msa, int N_seq, int L_msa);
int Find_prot(char *name, struct Prot_input *Prot, int *index, int N);
float **Read_function(char *file_fun, struct Prot_input *Prot,
		      int *index, int N);
void Sec_str_AA_propensity(struct protein *pi, struct protein *pj,
			   int *ali, int nali, float DP);
void Print_propensities(char *nameout);
double ***Prop_secstr=NULL, ***Prop_AA=NULL, *DP_bin=NULL, *DP_norm=NULL;
char *SS_name[4], AA_name[22];
extern int Code_AA_2(char aseq, char *aacode, int n);

// Auxiliary
void Initialize_prot(struct protein *prot);
float Normalization(float *norm_c,struct protein *proti,struct protein *protj);
void Write_ali_pair(int ***Ali_pair, int i, int j, int *L_seq, int *ali_PC);
char **Assign_Seq(char ***name_seq, int *len_seq,
		  struct Prot_input *Prot_in, int N_seq,
		  int *rep_str, int *i_seq, int N_ali);
int ***Select_alis(int ***Ali_pair, int N_seq, int *len_seq,int *rep_str);
void Set_contact_type();
int Count_AA(char *seq, int N_ali);
int *Representative_structure(int *N_conf, int **conformation, int n,
			      struct protein *prots, float **PC_all);
int Seq_differences(int *id, char *seq1, char *seq2, int N_ali);
int Count_gaps(char *seq1, char *seq2, int N_ali);
float Min_dist(int i, int j, int **conformation, int *N_conf, float **div);
float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div);
//void Change_conformations(int **conformation,int N_seq,int *N_conf,int *ali_str);
void Get_file_name(char *name, char *file);
float Divergence(float SI, float S0, float DMAX);
void Print_seq(int *ali, char *seq, int N_ali);
void Write_identity(int type,
		    struct protein *proti, struct protein *protj,
		    int *ali_ij, int *id_aa, int *id_sup, int *shift,
		    int *neigh_ali,int *neigh_noali,
		    int *neigh_noali_aaid,
		    int *ali_cont_ct, int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid, int *id_cont_aaid);
void Summary_identical(FILE *file_id, int it_opt);
void Count_noali(int *noali, int *ali,  int *ncont, int n, float c_ave);
float Mean_freq(double *se, double *re2, double sum, double tot);
int Count_ali(int *ali, int len);
double Sum_bins(double *bin);
float Ave_se(double *se, double sum1, double sum2, int n, float nind);
FILE *Open_file_div(char *name_div, char *name_in, char *head);
FILE *Open_file_sim(char *name_sim, char *name_in, char *head,
		    float **fun_sim);
void Print_scores(char *out, double *s1, char *what);

// Eliminate pairs with SI<SI_thr=2*S0=0.12
float SI_thr=0.10; //0.12

void Print_ali_ss(int *ali_i, short *ss_i, int *ali_j, short *ss_j, int N_ali);
char PDBDIR2[100]="";

/*********************************************************************
                          MAIN routine
**********************************************************************/
int main(int argc, char **argv)
{
  double time0=clock();

  // INPUT
  NORM=0; // Normalize with minimum (0) maximum (1) or geometric mean (2)?
  NORMA=0; // Normalize with nali (1) or with NORM (0)?
  int ITMAX1=ITMAX-1;
  char PDBDIR[100]="", PDB_EXT[10]=".pdb", OUTG[10]="PC";
  char file_ali[200]="", file_list[200]="", file_fun[200]="", name_in[80]="";
  Get_input(file_ali, file_list, file_fun, name_in, PDBDIR, PDB_EXT,
	    OUTG, &NORM, &ALI_SS, &SHIFT_MAX, &SS_MULT,
	    &PRINT_SIM, &PRINT_CV, &PRINT_DIV,
	    //&PRINT_PAIR, &PRINT_CLIQUE,
	    argc, argv);
  if(INP_MSA && pdblist[0]!='\0'){
    printf("ERROR, you cannot input both a list of pdb (-pdblist) and a "
	   "multiple alignment (-ali), please modify it\n"); exit(8);
  }
  if(INP_MSA==0 && ALI_SS){
    printf("WARNING, the input alignment is not an MSA and "
	   "it cannot be modified\nSetting ALI_SS=0\n"); ALI_SS=0; 
  }
  
  float d02=0; 
  
  // Input sequences and structures
  int N_ali=0, i, j, N_prot=0, N_prot_seq=0, N_ali_seq=0;
  struct Prot_input *Prot_in=NULL, *Prot_seq=NULL;
  if(file_list[0]!='\0'){
    N_prot=Get_pdb_list(&Prot_in, file_list);
  }
  if(file_ali[0]!='\0'){
    if(N_prot==0){
      N_prot=Get_sequences(&Prot_in, &N_ali, file_ali, INP_MSA);
    }else{
      N_prot_seq=Get_sequences(&Prot_seq, &N_ali_seq, file_ali, 0);
    }
  }
  
  // Output alignment types
  for(i=0; i<NTYPE; i++)code[i]=-1;
  ATYPE=0;
  if(INP_MSA){Set_type(ali_name, code, &ATYPE, "INP", 0);} // Input ali;
  else{Set_type(ali_name, code, &ATYPE, "PAIR", 0);}
  printf("INP_MSA= %d\n", INP_MSA);
  if(ALI_SS)Set_type(ali_name, code, &ATYPE, "SS", 1);
  if(ALI_TM)Set_type(ali_name, code, &ATYPE, "TM", 2);
  if(ALI_CO)Set_type(ali_name, code, &ATYPE, "CO", 3);
  if(ALI_PC)Set_type(ali_name, code, &ATYPE, "PC", 4);
  PTYPE=ATYPE; // Only pairwise alignments
  PT1=PTYPE-1;
  printf("PTYPE= %d\n", PTYPE);
  if(ALI_PC){
    if(NTYPE < PTYPE+NMSA){
      printf("WARNING, the number of built MSA must be <= PTYPE+%d=%d "
	     "but it is %d\n", NMSA, PTYPE+NMSA, NTYPE);
      NMSA=NTYPE-PTYPE;
      printf("Reduce MSA to %d\n", NMSA);
    }
    for(i=0; i<=NMSA; i++){
      char name[15]; sprintf(name, "PC_M%d", i);
      Set_type(ali_name, code, &ATYPE, name, 6+i);
    }
  }
  if(ALL_TYPES){CTYPE=PTYPE+1;}
  else{CTYPE=2;}
  
  float CO[ATYPE], TM[ATYPE], SI[ATYPE], PC[ATYPE], nali[ATYPE];
  double CO_sum[ATYPE], TM_sum[ATYPE], SI_sum[ATYPE],
    PC_sum[ATYPE], nali_sum[ATYPE];
  for(i=0; i<ATYPE; i++){
    opt_PC[i]=0; Diff_opt[i]=0;
    CO[i]=0; TM[i]=0; SI[i]=0; nali[i]=0; PC[i]=0;
    CO_sum[i]=0; TM_sum[i]=0; SI_sum[i]=0; nali_sum[i]=0; PC_sum[i]=0;
  }
  
  // What to print ? 
  //if((PRINT_CV==0)&&(PRINT_SIM==0))PRINT_DIV=1;
  
  
  /**************   READ PROTEIN STRUCTURES  ******************/
  // Read PDB files and compute contact matrices
  Set_contact_type();
  printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	 CONT_TYPE, CONT_THR,IJ_MIN);
  
  // Add final / to PDBDIR
  if(PDBDIR[0]!='\0'){
    char *c=PDBDIR, *c1=c+1;
    while(*c1!='\0'){c=c1; c1++;}
    if(*c!='/')strcat(PDBDIR,"/");
    printf("Looking for PDB files in directory %s\n", PDBDIR);
  }
  if(PDBDIR2[0]!='\0'){
    char *c=PDBDIR2, *c1=c+1;
    while(*c1!='\0'){c=c1; c1++;}
    if(*c!='/')strcat(PDBDIR2,"/");
    printf("Looking for PDB files in directory %s\n", PDBDIR2);
  }
  
  // Proteins prots
  struct protein prots[N_prot], *prot=prots;
  int N_pdb=0, i_seq[N_prot], i_inp[N_prot], L_seq[N_prot];
  int **Prot_ali=NULL; if(N_ali){Prot_ali=Allocate_mat2_i(N_prot, N_ali);}
  int n_no=0, n_nuc=0;

  for(int ir=0; ir<N_prot; ir++){
    struct Prot_input *Pn=Prot_in+ir;
    printf("\nReading PDB %s chain=%c (%d) %s %s\n",
	   Pn->code, Pn->chain, Pn->chain_num, Pn->domain, Pn->domname);
    Initialize_prot(prot);
    char *PDB_D=PDBDIR; if(Pn->dir){PDB_D=PDBDIR2;}
    int r=Read_PDB_compress(&prot, Pn->code, &(Pn->chain), Pn->chain_num,
			    PDB_D, PDB_EXT);
    if(r==0){
      r=Read_PDB_compress(&prot, Pn->domname, &(Pn->chain), Pn->chain_num,
			  PDB_D, PDB_EXT);
    }
    if(PDBDIR2[0]!='\0' && strcmp(PDB_D, PDBDIR2)){
      if(r==0){
	r=Read_PDB_compress(&prot, Pn->code, &(Pn->chain), Pn->chain_num,
			    PDBDIR2, PDB_EXT);
      }
      if(r==0){
	r=Read_PDB_compress(&prot, Pn->code, &(Pn->chain), Pn->chain_num,
			    PDBDIR2, PDB_EXT);
      }
    }
    if(r<=0){n_no++; continue;}
    if(prot->nuc){n_nuc++; Empty_prot(prot); continue;}
    
    strcpy(prot->domname, Pn->domname);
    strcpy(prot->domain,  Pn->domain);
    if(Pn->domain[0]!='\0'){
      Pn->nfrag=Read_domain(&(Pn->ini_frag), &(Pn->end_frag), Pn->domain);
      r=Select_domain(prot, Pn->ini_frag, Pn->end_frag, Pn->nfrag);
      if(r<=0)continue;
    }
    printf("L= %d\n", prot->len);
    
    if(Prot_in[ir].seq){ // Amino acids from alignment or FASTA
      Prot_in[ir].len=Count_AA(Prot_in[ir].seq, N_ali);
      if(Prot_in[ir].len!=prot->len){
	printf("WARNING, different n.residues in ali (%d) and PDB (%d)\n",
	       Prot_in[ir].len, prot->len);
	for(j=0; j<N_ali; j++){
	  if(Prot_in[ir].seq[j]!='-')printf("%c",Prot_in[ir].seq[j]);
	}
	printf("\n");
	for(j=0; j<prot->len; j++){printf("%c",prot->aseq[j]);}
	printf("\n");
      }
      if(Prot_ali && Align_seq(Prot_ali[N_pdb],N_ali,Prot_in[ir].seq,
			       prot->aseq,prot->len)<0)continue;
      if(Prot_in[ir].len<prot->len){prot->len=Prot_in[ir].len;}
    }else{ // Amino acids from PDB
      Prot_in[ir].len=prot->len;
      Prot_in[ir].seq = malloc(prot->len*sizeof(char));
      for(j=0; j<prot->len; j++)Prot_in[ir].seq[j]=prot->aseq[j];
    }
    
    L_seq[N_pdb]=prot->len;
    if(Prot_ali==NULL && prot->len>N_ali){N_ali=prot->len;}
    
    int NC=Compute_contact_list(prot, CONT_TYPE, CONT_THR, IJ_MIN);
    printf("%d contacts\n", NC);
    Set_scores(prot);
    i_seq[N_pdb]=N_pdb; i_inp[N_pdb]=ir; N_pdb++; prot++;
  }
  printf("\n%d proteins read out of %d listed in %s\n",N_pdb,N_prot,file_list);
  printf("%d chains not found, %d are nucleic acids total: %d\n",
	 n_no, n_nuc, n_no+n_nuc+N_pdb);

  if(N_pdb<2){
    printf("ERROR, fewer than 2 proteins found\n"); exit(8);
  }
  // End reading PDB files
  
  /*************************************************************************
       Group conformations of the same sequence
  **************************************************************************/
  // Group identical sequences
  // Structural divergence: minimum among all conformations, computed after
  // performing the multiple alignment only for non-redundant conformations 
  printf("Grouping proteins with > %.3f seq.identity\n", seq_thr);
  float str_thr=0.97; //**conformation, 
  int *rep_str, *N_all, **rep_conf, **rep_index, *N_conf,
    ***all_conf, **N_conf_all;
  int ***ali_clus, *L_ali_seq, ***ali_pairs;
  float ***TM_pairs;
  int N_seq=Cluster_seqs(&rep_str, &N_all, &rep_conf, &rep_index, &N_conf,
			 &all_conf, &N_conf_all,
			 &ali_clus, &L_ali_seq, &ali_pairs, &TM_pairs,
			 prots, N_pdb, seq_thr, str_thr);
  // rep_str[i]: representative structure of sequence i
  // N_conf[i] number of different conformations of sequence i
  //rep_conf[i][j]: representative conformation j of sequence i
  int num_conf=0; for(i=0; i<N_seq; i++){num_conf+=N_conf[i];}
  printf("%d structures grouped into %d sequences and %d conformations\n",
	 N_pdb, N_seq, num_conf);
  printf("Average conf. per sequence: %.1f\n", (float)num_conf/N_seq);
  printf("Average str. per conformation: %.1f\n", (float)N_pdb/num_conf);
 
  if(N_seq<=1){
    printf("ERROR, only one sequence found, rerun with higher seq_thr\n");
    exit(8);
  }
  
  int is_rep[N_pdb], is_conf[N_pdb], nconf=0;
  for(i=0; i<N_pdb; i++){is_rep[i]=-1; is_conf[i]=-1;}
  for(i=0; i<N_seq; i++){
    if(rep_str[i]>=0){is_rep[rep_str[i]]=i;}
    for(j=0; j<N_conf[i]; j++){
      if(rep_conf[i][j]>=0){is_conf[rep_conf[i][j]]=i; nconf++;}
    }
  }
  

  
  if(file_list[0]!='\0'){

    char name_new[200]; FILE *file_out;
    sprintf(name_new, "%s.nr.pdblist", name_in);
    file_out = fopen(name_new, "w");
    printf("Writing non-redundant seq in %s\n", name_new);
    for(i=0; i<N_pdb; i++){
      if(is_rep[i]<0){continue;}
      struct Prot_input *Pn=Prot_in+i_inp[i];
      fprintf(file_out, "%s %c %s %s 1\n",
	      Pn->code, Pn->chain, Pn->domain, Pn->domname);
    }
    fclose(file_out);
    
    if(nconf>N_seq){
      sprintf(name_new, "%s.nrconf.pdblist", name_in);
      file_out = fopen(name_new, "w");
      printf("Writing non-redundant conformations in %s\n", name_new);
      fprintf(file_out, "#pdb chain domain domname folder sequence\n");
      for(i=0; i<N_pdb; i++){
	if(is_conf[i]<0){continue;}
	struct Prot_input *Pn=Prot_in+i_inp[i];
	fprintf(file_out, "%s %c %s %s %d %d\n",
		Pn->code, Pn->chain, Pn->domain, Pn->domname, Pn->dir+1,
		is_conf[i]);
      }
      fclose(file_out);
    }
    if(PRINT_FASTA){

      // Write read sequences from PDB to FASTA file
      char name_fasta[200];
      sprintf(name_fasta, "%s.input.fasta", name_in);
      FILE *fasta_out = fopen(name_fasta, "w");
      printf("Writing input sequences to FASTA %s\n", name_fasta);
      for (int i = 0; i < N_pdb; i++) {
	struct protein *prot = prots + i;
	fprintf(fasta_out, ">%s\n", prot->domname);
	for (int j = 0; j < prot->len; j++) {
	  fprintf(fasta_out, "%c", prot->aseq[j]);
	}
	fprintf(fasta_out, "\n");
      }
      fclose(fasta_out);
    }
    if(EXIT_FASTA){
      printf("-exit_fasta was requested, exiting\n");
      exit(8);
    }
  }


  double time1=clock();
  printf("Identical sequences grouped.Time= %.2lf sec.\n",
	 (time1-time0)/nbtops);
 
  // Average number of contacts per residue
  c_ave=0; double c_norm=0;
  for(i=0; i<N_pdb; i++){
    struct protein *proti=prots+i;
    for(j=0; j<proti->len; j++)c_ave+=proti->ncont[j];
    c_norm+=proti->len;
  }
  c_ave/=c_norm;
  
  // Alignments
  if(Prot_ali==NULL){N_ali*=2; Prot_ali=Allocate_mat2_i(N_seq, N_ali);}
  int N_ali_max=3*N_ali, N_ali_max_ini=N_ali_max;
  
  /**************************************
      Secondary structure based alignment
  *******************************************/
  int **Prot_ali_ss=NULL, **Prot_ali_ij=NULL, NS=0;
  struct protein **prot_ij=NULL;
  
  if(Set_sec_str(prots, N_pdb)<0){ // Change - into c
    printf("WARNING, secondary structure information not found\n");
    ALI_SS=0; goto end_ss;
  }
  prot_ij=malloc(N_seq*sizeof(struct protein *));
  for(i=0; i<N_seq; i++)prot_ij[i]=prots+rep_str[i];
  //Write_ss_ali(Prot_ali,prot_ij,N_seq,N_ali,name_in, "Input");

  if(ALI_SS){
    
    int s_not=Test_notali(Prot_ali,prots,N_seq,N_ali);
    if(s_not)printf("WARNING, %d proteins had no aligned residues\n",s_not);
    if(SS_MULT){NS=N_seq;}else{NS=2;}
    Prot_ali_ss=malloc(NS*sizeof(int *));
    for(i=0; i<NS; i++)Prot_ali_ss[i]=malloc(N_ali*sizeof(int));
    
    if(SS_MULT){
      int N_ali_ss=
	Align_ss(Prot_ali_ss, Prot_ali, N_ali, prot_ij, N_seq, SHIFT_MAX, 1);
      Write_ss_ali(Prot_ali_ss, prot_ij,N_seq,N_ali_ss,name_in, "SSAli");
      Write_ali_prot(Prot_ali_ss,prot_ij, N_seq,N_ali_ss,name_in,"SSAli");
    }else{
      free(prot_ij);
      prot_ij=malloc(2*sizeof(struct protein *));
      Prot_ali_ij=malloc(2*sizeof(int *));
    }
  }else{ //ALI_SS==0
    free(prot_ij); prot_ij=NULL;
  }
  
  // Prepare output
 end_ss:
  printf("ALI_SS= %d SS_MULT= %d NS= %d SHIFT_MAX=%d\n",
	 ALI_SS,SS_MULT,NS,SHIFT_MAX);
  
  char ss_def[100];
  if(ALI_SS){
    sprintf(ss_def, "MSA corrected for secondary structure.");
    if(SS_MULT){strcat(ss_def, " Multiple alignment correction.");}
    else{strcat(ss_def, " Pairwise correction.");}
    char tmp[30]; sprintf(tmp, " SHIFT_MAX= %d\n",SHIFT_MAX);
    strcat(ss_def, tmp);
  }
  char norm_def[200];
  sprintf(norm_def,"# Normalization of nali: ");
  if(NORM==0){strcat(norm_def, " Minimum length\n");}
  else if(NORM==1){strcat(norm_def, " Maximum length\n");}
  else{strcat(norm_def, " Geometric mean\n");}
  if(NORMA){
    strcat(norm_def,"# Normalization of seqid, TM-score and CO: ");
    strcat(norm_def, " aligned residues\n");
  }
  
  char name_sim[100]; FILE *file_sim=NULL;
  if(PRINT_SIM_ALL){
    Change_ext(name_sim, name_in, EXT_SIM);
    file_sim=fopen(name_sim, "w");
    fprintf(file_sim, "### Prot1 Prot2 align Seq_Id Cont_Overlap TM_Score");
    if(ALL_TYPES && ALI_PC)
      fprintf(file_sim, " align_PC Seq_Id_PC Cont_Ov_PC TM_PC\n");
    fprintf(file_sim, "### 0 0  1 1 1 1");
    if(ALL_TYPES && ALI_PC)fprintf(file_sim, "  1 1 1 1");
    fprintf(file_sim, "\n");
  }
  
  // Allocate pairwise computations only for i>j
  
  ///////////////////////////////////
  int *ali_ij=malloc(N_ali_max*sizeof(int));  // Input MSA
  int *ali_tmp=malloc(N_ali_max*sizeof(int));
  int al2i=-1, al2j=-1;
  
  float **Seq_diff=malloc(N_seq*sizeof(float *));
  for(i=0; i<N_seq; i++){ // N_pdb
    Seq_diff[i]=malloc(N_seq*sizeof(float));
    Seq_diff[i][i]=0;
    L_seq[i]=prots[rep_str[i]].len;
    for(j=0; j<i; j++)Seq_diff[i][j]=L_seq[i];
  }
  
  int *ali_all[ATYPE];
  for(int it=0; it<ATYPE; it++){
    ali_all[it]=malloc(N_ali_max*sizeof(int));
  }
  
  // For computing divergence, PC is computed without tha nali part
  int RTYPE=3; // Reduced number of types. 0=Input 1=Optimal 2=TMP
  float **na_all[RTYPE], **SI_all[RTYPE], **TM_all[RTYPE],
    **CO_all[RTYPE], **PC_all[RTYPE];
  for(int it=0; it<RTYPE; it++){
    na_all[it]=malloc(N_seq*sizeof(float *));
    SI_all[it]=malloc(N_seq*sizeof(float *));
    TM_all[it]=malloc(N_seq*sizeof(float *));
    CO_all[it]=malloc(N_seq*sizeof(float *));
    PC_all[it]=malloc(N_seq*sizeof(float *));
    for(int i=0; i<N_seq; i++){
      int k; if(i){k=i;}else{k=1;}
      na_all[it][i]=malloc(k*sizeof(float));
      SI_all[it][i]=malloc(k*sizeof(float));
      TM_all[it][i]=malloc(k*sizeof(float));
      CO_all[it][i]=malloc(k*sizeof(float));
      PC_all[it][i]=malloc(k*sizeof(float));
    }
  }

  // Conditional probabilities
  if(PRINT_ID){
    for(i=0; i<ATYPE; i++){
      id_aa[i]=malloc(N_ali*sizeof(int));
      id_sup[i]=malloc(N_ali*sizeof(int));
      neigh_ali[i]=malloc(N_ali*sizeof(int));
      neigh_noali[i]=malloc(N_ali*sizeof(int));
      neigh_noali_aaid[i]=malloc(N_ali*sizeof(int));
      shift[i]=malloc(N_ali*sizeof(int));
      id_cont_ct[i]=malloc(2*sizeof(int));
      ali_cont_ct[i]=malloc(2*sizeof(int));
      id_cont_sup[i]=malloc(3*sizeof(int));
      ali_cont_sup[i]=malloc(3*sizeof(int));
      id_cont_aaid[i]=malloc(3*sizeof(int));
      ali_cont_aaid[i]=malloc(3*sizeof(int));
    }
  }

  // PC0
  PC0_1=PC_load[1]*S0+PC_load[2]*TM0; // excluding nali
  PC_norm_div=PC_load[1]+PC_load[2]+PC_load[3];
  PC_norm=PC_norm_div+PC_load[0];
  if(DIV_ALIGNED){PC0_1+=PC_load[0]*0.5; PC_norm_div+=PC_load[0];}
  // PC value for random pairs
  
  // Store PCA alignments
  //int **PC_opt_pair=Allocate_mat2_i(N_pdb, N_pdb);
  float **d2=Allocate_mat2_f(N_ali, N_ali);
  int   **nc=Allocate_mat2_i(N_ali, N_ali);
  int **Ali_pair[N_seq]; // ali[i][j][site_i]=site_j i<j
  float ***Rot_pair[N_seq], **Shift_pair[N_seq];
  for(i=0; i<N_seq; i++){ // Allocate for all pairs
    int n; if(ALL_PAIRS){n=N_seq;}else{n=i; if(n==0)continue;}
    Ali_pair[i]=Allocate_mat2_i(n, L_seq[i]);
    Shift_pair[i]=Allocate_mat2_f(n, 3);
    Rot_pair[i]=malloc(n*sizeof(float **));
    for(j=0; j<n; j++)Rot_pair[i][j]=Allocate_mat2_f(3,3);
  }
  
  
  /*************************************************************
         Multiple superimposition
  *************************************************************/
  struct protein *prot_p[N_seq]; 
  for(i=0; i<N_seq; i++){prot_p[i]=prots+rep_str[i];}
  float TMs=-1;
  if(INP_MSA){    // Multiple structure alignment
    TMs=TM_score_mult(TM_all[0], Prot_ali, N_ali, N_seq, prot_p, 0);
    if(TMs>0){
      printf("Multiple superimposition done, TM= %.3f\n",
	     TMs*2/(N_seq*(N_seq-1)));
    }else{
      printf("WARNING, multiple superimposition could not be done\n");
    }
  }
  
  /*************************************************************
         Pairwise computations for all PDB pairs j<i
  *************************************************************/
  printf("Computing pairwise structural scores\n");
  int seq_group[N_seq];
  for(i=0; i<N_seq; i++){
    int al1i=i_seq[i];
    struct protein *proti=prots+rep_str[i];
    printf("Seq %s %d of %d L= %d\n", proti->domname, i, N_seq, proti->len);
    printf("Pairwise alignments: ");
    seq_group[i]=-1;
    
    if(ALI_SS && SS_MULT==0){
      prot_ij[0]=prots+i; Prot_ali_ij[0]=Prot_ali[i];
    }
    
    for(j=0; j<i; j++){
      struct protein *protj=prots+rep_str[j];
      npair_pdb++;
      
      // Sequence differences
      int id, al1j=i_seq[j];
      
      // Normalization
      float norm_c, norm_ali=Normalization(&norm_c, proti, protj);
      
      // General
      int Comp_TM=1;
      
      // Input alignment
      int it=0;
      if(INP_MSA){
	Pair_ali(ali_ij, N_ali, Prot_ali[i], Prot_ali[j])/norm_ali;
	if(TMs>0){
	  TM[0]=TM_all[0][i][j];
	  //Comp_TM=0;
	}
      }else{ // Perform initial pairwise alignment
	//printf("Performing pairwise alignment\n");
	//Align_PC_NW_pair(ali_ij, NULL, NULL, proti, protj);
	Invert_ali(ali_ij, proti->len, ali_pairs[i][j], protj->len);
	//int *tmp=ali_pairs[i][j];  // i1=i>i2=j
	//for(int k=0; k<proti->len; k++){ali_ij[k]=tmp[k];}
      }
      
      Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff, 
		      na_all[0],SI_all[0],TM_all[0],CO_all[0],PC_all[0],
		      NULL, NULL, c_ave, &d02, d2, ali_all[it],
		      ali_ij, Comp_TM, TMs, proti, protj,
		      norm_c, norm_ali, it, 0);
      Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		 nali, SI, TM, CO, PC, it);

      /*// Leave if sequence identity is above threshold
	if(OMIT_SAME && SI[it] > seq_thr){
	printf("prots %s and %s are almost identical, omitting %s\n",
	proti->domname, protj->domname, proti->domname);
	seq_group[i]=seq_group[j];
	goto end_i;
	}*/
      
      Comp_TM=1;
      // Sec.str. corrected alignments 
      if(ALI_SS){
	it++;
	// Pairwise alignment
	if(SS_MULT==0){
	  prot_ij[1]=protj; Prot_ali_ij[1]=Prot_ali[j];
	  Align_ss(Prot_ali_ss,Prot_ali_ij,N_ali,prot_ij,NS,SHIFT_MAX,0);
	  Pair_ali(ali_ij, N_ali, Prot_ali_ss[0], Prot_ali_ss[1]);
	}else{
	  Pair_ali(ali_ij, N_ali, Prot_ali_ss[i], Prot_ali_ss[j]);
	}
	Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			na_all[2],SI_all[2],TM_all[2],CO_all[2],PC_all[2],
			NULL, NULL, c_ave, &d02, d2, ali_all[it],
			ali_ij, Comp_TM, TMs, proti, protj,
			norm_c, norm_ali, it, 0);
      Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		 nali, SI, TM, CO, PC, it);

	int accept_SS=0, SI_high=-1, TM_high=-1, CO_high=-1;
	if(SI[it]<SI[0]){SI[it]=SI[0]; SI_high=0;}
	else if(SI[it]>SI[0]){SI_high=1; }
	if(SI_high>=0)sum_SS_SI[SI_high]++;
	if(TM[it]>=TM[0]){
	  accept_SS=1; if(TM[it]>TM[0]){TM_high=1;}
	}else{
	  TM[it]=TM[0]; TM_high=0;
	}
	if(TM_high>=0)sum_SS_TM[TM_high]++;

	if(0 && (TM[it]-TM[0])<-0.03){
	  int l=proti->len; if(protj->len<l)l=protj->len;
	  printf("WARNING, Decrease of TM from %.1f to %.1f\n",
		 TM[0]*l,TM[it]*l);
	  Print_ali_ss(Prot_ali[i], proti->ss3,
		       Prot_ali[j], protj->ss3, N_ali);
	  printf("=========================================\n");
	  int ii=0,jj=1; if(SS_MULT){ii=i; jj=j;}
	  Print_ali_ss(Prot_ali_ss[ii], proti->ss3,
		       Prot_ali_ss[jj], protj->ss3, N_ali);
	}

	if(CO[it]<CO[0]){CO_high=0;}
	else if(CO[it]>CO[0]){CO_high=1;}
	if(accept_SS==0)CO[it]=CO[0];
	if(CO_high>=0)sum_SS_CO[CO_high]++;
	if(SI_high>=0 && TM_high>=0)sum_SS_SI_TM[SI_high][TM_high]++;
	if(SI_high>=0 && CO_high>=0)sum_SS_SI_CO[SI_high][CO_high]++;
	if(TM_high>=0 && CO_high>=0)sum_SS_TM_CO[TM_high][CO_high]++;
      }

      // TM_score based alignment
      float d2min1[proti->len];
      if(ALI_TM){
	it++; int *ali_TM=ali_all[it], *ali_ini=ali_all[it-1], k;
	for(k=0; k<proti->len; k++)ali_tmp[k]=ali_ini[k];
	for(int iter=0; iter<ITMAX; iter++){
	  Align_TM(ali_TM, d2min1, d2, d02, ali_tmp, proti->len, protj->len);
	  TM[it]=TM_score(d2, &d02, NULL, NULL, ali_TM, norm_ali,
	  		  proti, protj, 0);
	  for(k=0; k<proti->len; k++)ali_tmp[k]=ali_TM[k];
	}
	Comp_TM=0;
	Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			na_all[2],SI_all[2],TM_all[2],CO_all[2],PC_all[2],
			NULL, NULL, c_ave, &d02, d2, ali_all[it],
			ali_TM, Comp_TM, TMs, proti, protj,
			norm_c, norm_ali, it, 0);
	Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		   nali, SI, TM, CO, PC, it);
      }

     // Contact Overlap based alignment
      if(ALI_CO){
	it++; int *ali_CO=ali_all[it], *ali_ini=ali_all[it-1], k;
	for(k=0; k<proti->len; k++)ali_tmp[k]=ali_ini[k];
	for(int iter=0; iter<ITMAX; iter++){
	  Align_CO(ali_CO, ali_tmp, nc,
		   proti->Cont_map, proti->len, protj->Cont_map, protj->len);
	  CO[it]=
	    Contact_overlap(ali_CO, 
			    proti->Cont_map, proti->len,
			    protj->Cont_map, protj->len);
	  for(k=0; k<proti->len; k++)ali_tmp[k]=ali_CO[k];
	}

	Comp_TM=1; 
	Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			na_all[2],SI_all[2],TM_all[2],CO_all[2],PC_all[2],
			NULL, NULL, c_ave, &d02, d2, ali_all[it],
			ali_CO, Comp_TM, TMs, proti, protj,
			norm_c, norm_ali, it, 0);
	Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		   nali, SI, TM, CO, PC, it);
      }


     // PC based alignment
      if(ALI_PC){
	it++; int *ali_PC=ali_all[it];
	int *ali_ini=ali_all[it-1], k;

	for(int iter=0; iter<ITMAX; iter++){
	  if(PC_NW){
	    Align_PC_NW_pair(ali_PC, d2, nc, proti, protj);
	  }else{
	    if(iter==0)for(k=0; k<proti->len; k++)ali_tmp[k]=ali_ini[k];
	    Align_PC(ali_PC, ali_tmp, PC_load, d02, d2, nc,
		     proti->aseq, proti->len, protj->aseq, protj->len);
	    for(k=0; k<proti->len; k++)ali_tmp[k]=ali_PC[k];
	  }
	  if(iter==ITMAX1){
	    TM[it] = TM_score(d2, &d02, Rot_pair[i][j], Shift_pair[i][j],
			      ali_PC, norm_ali, proti, protj, 0);
	    if(0)Test_Rot(d2,i,j,
			  Rot_pair[i][j],Shift_pair[i][j],
			  proti->xca_rot, proti->len,
			  protj->xca_rot, protj->len);
	  }else{
	    TM[it] = TM_score(d2, &d02, NULL, NULL,
			      ali_PC, norm_ali, proti,protj,0);
	    Shared_contacts(nc, ali_PC,
			    proti->Cont_map, proti->len,
			    protj->Cont_map, protj->len);
	  }
	}
	
	Comp_TM=PRINT_ID; 
	Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			na_all[2],SI_all[2],TM_all[2],CO_all[2],PC_all[2],
			NULL, NULL, c_ave, &d02, d2, ali_all[it],
			ali_PC, Comp_TM, TMs, proti, protj,
			norm_c, norm_ali, it, 0);
	Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		   nali, SI, TM, CO, PC, it);
      }
      
      // Optimal alignments
      int k_PC=-1;
      for(it=0; it<PTYPE; it++){
	int ko=0, k; // Optimal alignment
	if(ALI_PC && strcmp(ali_name[it],"PC")==0){ // PC_ali
	  int kk=0; 
	  for(k=1; k<PTYPE; k++)if(k!=it && PC[k]>PC[kk])kk=k;
	  opt_PC[kk]++;
	  if(PC[it]>=PC[kk]){ko=it;}else{ko=kk;}
	  Write_ali_pair(Ali_pair, i, j, L_seq, ali_all[ko]);
	  //printf("Writing alignment %d %d\n", i, j);
	  printf(" %d", j);
	  k_PC=ko;
	  //PC_opt_pair[i][j]=ko;
	}else if(ALI_CO && strcmp(ali_name[it],"CO")==0){ // CO_ali
	  for(k=1; k<PTYPE; k++)if(CO[k]>CO[ko])ko=k;
	}else if(ALI_TM && strcmp(ali_name[it],"TM")==0){ // TM_ali
	  for(k=1; k<PTYPE; k++)if(TM[k]>TM[ko])ko=k;
	}else{
	  continue;
	}
	if(ko!=it){
	  Diff_opt[it]++;
	}
      }
      
      // Print similarities
      if(0 && file_sim){
	fprintf(file_sim, "%s\t%s", proti->domname, protj->domname);
	for(int a=0; a<PTYPE; a++){
	  if(a){
	    if(ALL_TYPES==0)break; // Print only input
	    if(ALL_TYPES && a< PT1)continue; // Print input+ALI_PC
	  }
	  fprintf(file_sim,"\t%.3f\t%.3f\t%.3f\t%.3f",
		  nali[a],SI[a],CO[a],TM[a]);
	}
	fprintf(file_sim,"\n");
      }
    } // end j
    printf("\n"); // Print pairwise alignment on screen
    seq_group[i]=i;
  end_i: continue;
  } // end pairs
  Empty_matrix_i(nc, N_ali);
  Empty_matrix_f(d2, N_ali);

  double time2=clock();
  printf("End of pairwise alignments. Time= %.2lf sec.\n",
	 (time2-time1)/nbtops);

  printf("End pairwise computations\n");
  if(file_sim){
    printf("Similarities written in file %s for all PDB pairs\n", name_sim);
    fclose(file_sim); 
  }

  
  /*************************************************************
         End pairwise computations for all PDB pairs j<i
  *************************************************************/
  
  /*************************************************************
          Optimal pairwise alignment it_opt_p
  *************************************************************/
  int it_opt_p=0;
  for(int k=1; k<PTYPE; k++){if(opt_PC[k]>opt_PC[it_opt_p])it_opt_p=k;}
  printf("Optimal pairwise alignment: %d %d pairs\n",
	 it_opt_p, opt_PC[it_opt_p]);
  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){
      na_all[1][i][j]=na_all[2][i][j];
      SI_all[1][i][j]=SI_all[2][i][j];
      TM_all[1][i][j]=TM_all[2][i][j];
      CO_all[1][i][j]=CO_all[2][i][j];
      PC_all[1][i][j]=PC_all[2][i][j];
    }
  }

  char head[1000]="";
  if(ALI_SS){strcat(head, "# "); strcat(head, ss_def);}
  strcat(head, norm_def); 
  
  /********************************************************************
        PC divergence between reference structures WITHOUT scoring nali
  *********************************************************************/

  float **PC_Div=Allocate_mat2_f(N_seq, N_seq), **PC_Div_opt=NULL;
  double qinf=0; int homo;
  int it=2; //it_opt_p; //it=0; // Only for input pairwise alignment
  for(i=0; i<N_seq; i++){
    struct protein *proti=prots+rep_str[i];
    for(j=0; j<i; j++){
      struct protein *protj=prots+rep_str[j];
      Compute_Dcont(&qinf,CO_all[it][i][j],
		    proti->len,protj->len,&homo,NORM);
      double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm_div;;
      double PC_a = PC_all[it][i][j];
      if(DIV_ALIGNED==0){PC_a=(PC_a*PC_norm-na_all[it][i][j])/PC_norm_div;}
      PC_Div[i][j]=Divergence(PC_a, PC0, DMAX);
    }
  }
  printf("End of pairwise computations\n");
  
  /******************************************************************/
  /*                        Multiple alignments                     */
  /******************************************************************/
  
  int it_opt=-1;
  int **msa_opt=NULL, L_msa_opt=0;
  char name_PC[90];
  
  printf("Making %d multiple alignments\n", NMSA);
  int **msa=Allocate_mat2_i(N_seq, N_ali_max);
  msa_opt=Allocate_mat2_i(N_seq, N_ali_max);
  PC_Div_opt=Allocate_mat2_f(N_seq, N_seq);
  float **Seq_Div=Allocate_mat2_f(N_seq, N_seq);

  float PC_opt=0;
  int len_seq[N_seq], L_msa=0;
  char name_msa[95];
  sprintf(name_PC, "%s.PCAli", name_in);
  sprintf(name_msa, "%s.fas", name_PC);
  
  // Sequences
  char **name_seq, **Seq=
    Assign_Seq(&name_seq, len_seq, Prot_in,N_seq,rep_str,i_inp,N_ali);
  
  if(Prot_ali_ss){
    for(i=0; i<NS; i++){free(Prot_ali_ss[i]);} free(Prot_ali_ss);
  }
  
  // Pointers to proteins
  for(i=0; i<N_seq; i++)prot_p[i]=prots+rep_str[i];
  
  NMSA_act=NMSA;  
  npair_seq=N_seq*(N_seq-1)/2;
  int it_down=0, it_down_max=3;
  int last=0, itt_last=NMSA-1, Comp_TM=0;
  for(int itt=0; itt<=NMSA; itt++){
    int it=PTYPE+itt; 
    if(itt==0 && MAKE_CLIQUE){
      int N_ali_ini=N_ali; if(INP_MSA==0){N_ali_ini*=1.5;}
      Clique_MSA(msa, &L_msa, N_ali_max, Ali_pair, len_seq,
		 N_seq, Seq, name_seq, name_msa, N_ali_ini);
      
    }else if(last){
      // Restore optimal values
      L_msa=L_msa_opt;
      for(i=0; i<N_seq; i++){ // Copy msa and PC_Div
	int *mso=msa_opt[i], *ms=msa[i], j;
	for(j=0; j<L_msa; j++){*ms=*mso; mso++; ms++;}
	float *Div_o=PC_Div_opt[i], *Div_s=PC_Div[i];
	for(j=0; j<i; j++){*Div_s=*Div_o; Div_o++; Div_s++;}
      }
      
    }else{ // Progressive alignment
      for(i=0; i<N_seq; i++) 
	for(j=0; j<i; j++)PC_Div[j][i]=PC_Div[i][j]; // Symmetrize
      
      if(itt==0 || TMs<0){
	L_msa=NJ_align(msa, &N_ali_max, Ali_pair, Rot_pair, Shift_pair,
		       prot_p, name_PC, PC_Div, NULL, N_seq, AV_LINK);
      }else{
	L_msa=NJ_align(msa, &N_ali_max, Ali_pair, NULL, NULL,
		       prot_p, name_PC, PC_Div, NULL, N_seq, AV_LINK);
      }
      L_msa=Remove_gap_cols(msa, N_seq, L_msa);
    }
    
    // Score alignment
    //float **d2=NULL; //Allocate_mat2_f(L_msa, L_msa);
    printf("Score multiple alignment %d\n", itt);
    TMs=TM_score_mult(TM_all[PTYPE], msa, L_msa, N_seq, prot_p, 0);
    // The rotated coordinates are recorded at prot_p[i]
    if(TMs>0){
      printf("Multiple superimposition done, TM= %.3f\n",
	     TMs*2/(N_seq*(N_seq-1)));
      Comp_TM=1; //=0;
    }else{
      printf("WARNING, multiple superimposition could not be done\n");
      Comp_TM=1;
    }
    if(last){
      //it=it_opt;
      d2=Allocate_mat2_f(L_msa, L_msa);
      if(PRINT_ID)Comp_TM=1;
      printf("Last MSA (%d), Comp_TM=%d\n", it, Comp_TM);
    }else{
      d2=NULL;
    }
    if(N_ali_max>N_ali_max_ini){
      free(ali_ij);
      ali_ij=malloc(N_ali_max*sizeof(int));
      N_ali_max_ini=N_ali_max;
    }
    
    // Score pairwise
    for(i=0; i<N_seq; i++){
      struct protein *proti=prot_p[i]; //prots+rep_str[i];
      printf("Str %s %d of %d L= %d\n", proti->domname,i,N_seq,proti->len);
      for(j=0; j<i; j++){
	struct protein *protj=prot_p[j]; //prots+rep_str[j];
	Pair_ali(ali_ij, L_msa, msa[i], msa[j]);
	
	float d02, norm_c, norm_ali=Normalization(&norm_c, proti, protj);
	if(TMs>0){TM[it]=TM_all[PTYPE][i][j];}
	Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			na_all[2],SI_all[2],TM_all[2],CO_all[2],PC_all[2],
			Rot_pair[i][j], Shift_pair[i][j],
			c_ave, &d02, d2, Ali_pair[i][j], ali_ij, 
			Comp_TM, TMs, proti, protj,
			norm_c, norm_ali, it, last);
	Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
		   nali, SI, TM, CO, PC, it);
	double qinf=0; int homo;
	Compute_Dcont(&qinf,CO[it],proti->len,protj->len,&homo,NORM);
	double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm_div;
	double PC_a=PC[it];
	if(DIV_ALIGNED==0){PC_a=(PC_a*PC_norm-nali[it])/PC_norm_div;}
	PC_Div[i][j]=Divergence(PC_a, PC0, DMAX);
	if(last){ // Propensity of sec.str.
	  Sec_str_AA_propensity(proti, protj,ali_ij,nali[it],PC_Div[i][j]);
	}
      }
    } // end pairs

    if(it_opt<0 || PC_sum[it] > PC_opt){
      PC_opt=PC_sum[it]; L_msa_opt=L_msa; it_opt=it; NMSA_act=itt+1;
      for(i=0; i<N_seq; i++){ // Copy msa and PC_div
	int *mso=msa_opt[i], *ms=msa[i], j;
	for(j=0; j<L_msa; j++){*mso=*ms; mso++; ms++;}
	float *Div_o=PC_Div_opt[i], *Div_s=PC_Div[i];
	for(j=0; j<i; j++){*Div_o=*Div_s; Div_o++; Div_s++;}
	for(j=0; j<i; j++){
	  na_all[1][i][j]=na_all[2][i][j];
	  SI_all[1][i][j]=SI_all[2][i][j];
	  TM_all[1][i][j]=TM_all[2][i][j];
	  CO_all[1][i][j]=CO_all[2][i][j];
	  PC_all[1][i][j]=PC_all[2][i][j];
	  Seq_Div[i][j]=Divergence(SI_all[1][i][j], S0, DMAX);
	}
      }
      it_down=0;
    }else{
      it_down++;
      if(it_down==it_down_max){itt_last=itt; NMSA_act=itt+1;} 
    }
    if(last){break;}
    else if(itt==itt_last){last=1;}
    else{
      for(int it2=it-1; it2>=PTYPE; it2--){
	if(fabs(PC_sum[it]-PC_sum[it2])<0.00001*npair_seq){
	  printf("Score has repeated, exiting\n");
	  NMSA_act=itt+1; last=1; break;
	}
      }
    }
    
  } // end iter for finding optimal MSA
  double time3=clock();
  printf("Multiple alignments: %.2lf sec. \n",(time3-time2)/nbtops);

  // Print optimal MSA
  printf("Initial MSA length: %d Optimal: %d\n", N_ali, L_msa_opt);
  Print_MSA(msa_opt, name_msa,N_seq,L_msa_opt, len_seq, Seq, name_seq);
  Write_ss_ali(msa_opt, prot_p,N_seq, L_msa_opt, name_in, "PCAli");
  if(PRINT_PDB){Print_pdb(name_in, prot_p, N_seq);}


  if(ALL_CONF){
    /***********************************************************************
        Max PC divergence over all conformations of the same protein
    ************************************************************************/
    
    printf("Computing maximum structure similarity "
	   "over all conformations of the same protein\n");
    strcat(head, "# PC similarity maximized over all different "
	   "conformations of the same protein and different alignments\n");
    int it=1;  // optimal alignment

    for(i=0; i<N_seq; i++){
      int ci=rep_str[i];
      struct protein *proti=prots+ci;

      for(j=0; j<i; j++){
	if(N_conf[i]==1 && N_conf[j]==1){continue;}
	// Max PC over conformations of the same protein
	float PC_max=PC_all[it][i][j], PC_min=PC_max;

	int cj=rep_str[j];
	struct protein *protj=prots+cj;
	float norm_c, norm_ali=Normalization(&norm_c, proti, protj);

	for(int ki=0; ki<N_conf[i]; ki++){
	  int cki=rep_conf[i][ki], kki=rep_index[i][ki];
	  for(int kj=0; kj<N_conf[j]; kj++){
	    int ckj=rep_conf[j][kj], kkj=rep_index[j][kj];
	    if(cki==ci && ckj==cj){continue;}
	    struct protein *prot_ki=prots+cki;
	    struct protein *prot_kj=prots+ckj;
	    Pair_ali_conf(ali_ij, L_msa, msa_opt[i], msa_opt[j],
			  ali_clus[i][kki], L_ali_seq[i], prot_ki->len,
			  ali_clus[j][kkj], L_ali_seq[j], prot_kj->len);
	    int Comp_TM=1, TMs=0;

	    Score_alignment(nali, SI, TM, CO, PC, i, j, Seq_diff,
			    NULL, NULL, NULL, NULL, NULL, NULL, NULL,
			    c_ave, &d02, d2, ali_all[it],
			    ali_ij, Comp_TM, TMs, prot_ki, prot_kj,
			    norm_c, norm_ali, it, 0);
	    /*Sum_scores(nali_sum, SI_sum, TM_sum, CO_sum, PC_sum,
	      nali, SI, TM, CO, PC, it); */
	    if(PC[it]>PC_max){
	      //protj_opt=protj;
	      PC_max=PC[it];
	      na_all[it][i][j]=nali[it];
	      SI_all[it][i][j]=SI[it];
	      TM_all[it][i][j]=TM[it];
	      CO_all[it][i][j]=CO[it];
	      PC_all[it][i][j]=PC[it];
	    }
	    if(PRINT_CV && PC[it]<PC_min){ // Upper diagonal
	      PC_min=PC[it];
	      na_all[it][i][j]=nali[it];
	      SI_all[it][i][j]=SI[it];
	      TM_all[it][i][j]=TM[it];
	      CO_all[it][i][j]=CO[it];
	      PC_all[it][i][j]=PC[it];
	    }
	  } // end loop on conformations
	} // end loop on conformations
	
	// Compute PC_div without naligned
	double qinf=0; int homo;
	Compute_Dcont(&qinf,CO_all[it][i][j],
		      proti->len,protj->len,&homo,NORM);
	double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm_div;
	double PC_a = PC_max;
	if(DIV_ALIGNED==0){PC_a=(PC_a*PC_norm-na_all[1][i][j])/PC_norm_div;}
	PC_Div_opt[i][j]=Divergence(PC_a, PC0, DMAX);
      } // end loop of sequences
    }// end loop of sequences
  } //end ALL_CONF
  
    /***********************************************************************
        PC divergence for building the tree
    ************************************************************************/
  for(i=0; i<N_seq; i++) {
    PC_Div_opt[i][i]=0;
    for(j=0; j<i; j++)PC_Div_opt[j][i]=PC_Div_opt[i][j]; // Symmetrize
  }
  
  // Print PC Divergence
  if(PRINT_DIV){
    char name_div[100];
    Change_ext(name_div, name_in, ".PC.div");
    printf("Printing PC divergence in %s\n", name_div);
    FILE *file_div=fopen(name_div, "w");
    fprintf(file_div, "%d\n", N_seq);
    for(i=0; i<N_seq; i++){
      fprintf(file_div, "%s", prots[rep_str[i]].domname);
      for(j=0; j<N_seq; j++)fprintf(file_div, " %.3f",PC_Div_opt[i][j]);
      fprintf(file_div, "\n");
    }
    fclose(file_div);
  }
  
  /*********************************************************
           Distance-based Trees
  ***********************************************************/
  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){Seq_Div[j][i]=Seq_Div[i][j];}
  }
  NJ_align(NULL, &N_ali_max, NULL, NULL, NULL,
	   prot_p, name_PC, PC_Div_opt, Seq_Div, N_seq, 0); // NJ tree
  NJ_align(NULL, &N_ali_max, NULL, NULL, NULL,
	   prot_p, name_PC, PC_Div_opt, Seq_Div, N_seq, 1); // AL tree
  
  /********************* Print summary **********************/
  {
    char name_out[100]; sprintf(name_out, "%s.summary.dat", name_PC);
    printf("Writing summary scores in %s\n", name_out);
    FILE *file_out=fopen(name_out, "w");
    char out[2000]="";
    for(int k=0; k<5; k++){
      Sim_ave[k]=malloc(ATYPE*sizeof(double));
      Sim_diff[k]=Allocate_mat2_d(PTYPE+1, PTYPE+1);
      Sim_diff_2[k]=Allocate_mat2_d(PTYPE+1, PTYPE+1);
    }
    for(i=0; i<ATYPE; i++){
      Sim_ave[0][i]=nali_sum[i];
      Sim_ave[1][i]=SI_sum[i];
      Sim_ave[2][i]=TM_sum[i];
      Sim_ave[3][i]=CO_sum[i];
      Sim_ave[4][i]=PC_sum[i];
    }

    Print_scores(out, Sim_ave[0], "align length");
    Print_scores(out, Sim_ave[1], "seq identity");
    Print_scores(out, Sim_ave[2], "TM score");
    Print_scores(out, Sim_ave[3], "contact overlap");
    Print_scores(out, Sim_ave[4], "Princ Component");
    fprintf(file_out, "%s", out);
    fclose(file_out);
    printf("Optimal multiple alignment out of %d: %s %.4f\n",
	   NMSA, ali_name[it_opt], Sim_ave[4][it_opt]/npair_seq);
  }

  printf("MSA of %d PDB and %d different sequences\n", N_pdb, N_seq); 
  printf("Total time= %.2lf sec.\n",(clock()-time0)/nbtops);
  printf("Initial clustering: %.2lf sec.\n",(time1-time0)/nbtops);
  printf("Pairwise alignments: %.2lf sec.\n",(time2-time1)/nbtops);
  printf("Multiple alignments: %.2lf sec.\n", (time3-time2)/nbtops);
  printf("Conformations and printing: %.2lf sec.\n\n",(clock()-time3)/nbtops);
  
  /*****************************************************************
                 Print sec.str. propensities for multiple PC_ali
  ******************************************************************/
  printf("Printing propensities\n");
  char name_prop[100];
  sprintf(name_prop, "%s.id", name_PC);
  Print_propensities(name_prop);
  printf("Done\n");
  
  /****************************************************************
              Function similarity for different sequences
  *****************************************************************/
  // Functional similarity (if file_fun is present)
  float **fun_sim_Seq=NULL;
  int n_low=0, n_high=0;
  float **fun_sim=Read_function(file_fun, Prot_in, i_inp, N_pdb);
  if(fun_sim){
    fun_sim_Seq=Allocate_mat2_f(N_seq, N_seq);
    
    for(i=0; i<N_seq; i++){
      for(j=0; j<N_seq; j++){
	int ci=rep_str[i], cj=rep_str[j], cc;
	if(cj>ci){cc=ci; ci=cj; cj=cc;}
	if(fun_sim[ci][cj]<0)continue;
	float s=fun_sim[ci][cj];
	fun_sim_Seq[i][j]=s;
	fun_sim_Seq[j][i]=s;
	if(s<=fun_low)n_low++;
	if(s>=fun_high)n_high++;
      }
    }
    
    Empty_matrix_f(fun_sim, N_pdb);
    printf("%d pairs with function similarity <= %.3f\n",
	   n_low, fun_low);
    printf("%d pairs with function similarity >= %.3f\n",
	   n_high, fun_high);
  }
  
  // Optimal alignment
  int o=-1; // Scores based on PC_ALI, o=5
  if(it_opt>=0){o=it_opt;}
  else if(ALI_PC){o=code[4];}
  else if(ALI_CO){o=code[3];}
  else if(ALI_SS){o=code[1];}
  else if(ALI_TM){o=code[2];}
  if(o<0){
    printf("ERROR, optimal alignment undefined\n"); exit(8);
  }
  
  
  /****************************************************************
              Print pairwise similarity and divergence scores
  *****************************************************************/
  
  // Open output files
  char name_div_all[100];
  FILE *file_div_all=NULL; file_sim=NULL;
  //if(PRINT_DIV_ALL)file_div_all=Open_file_div(name_div_all, name_in, head);
  if(PRINT_SIM){file_sim=Open_file_sim(name_sim, name_in, head, fun_sim_Seq);}
  float **TN_Div=Allocate_mat2_f(N_seq, N_seq);
  float **TM_Div=Allocate_mat2_f(N_seq, N_seq);
  float **Cont_Div=Allocate_mat2_f(N_seq, N_seq);
  float **Seq_Id=Allocate_mat2_f(N_seq, N_seq);
  
  for(i=0; i<N_seq; i++){
    int ci=rep_str[i]; struct protein *pi=prots+ci;
    for(j=0; j<i; j++){
      int cj=rep_str[j]; struct protein *pj=prots+cj;
      
      if(file_sim)fprintf(file_sim, "%s\t%s", pi->domname, pj->domname);
      //if(file_div_all)fprintf(file_div_all, "%s\t%s",pi->domname,pj->domname);
      
      for(int it=0; it<2; it++){//it=0: input MSA it=1: opt MSA
	int a; 
	
	if(it==0){if(INP_MSA){a=0;}else{continue;}}
	else{
	  //if(ALL_TYPES){a=PTYPE;}else{a=1;}
	  a=1; //it_opt;
	}
	if(file_sim){
	  fprintf(file_sim, "\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f",
		  na_all[a][i][j], SI_all[a][i][j],CO_all[a][i][j],
		  TM_all[a][i][j], PC_all[a][i][j]);
	  if(msa_opt){
	    float RMSD=Compute_RMSD(pi, pj, i, j, msa_opt, L_msa_opt);
	    fprintf(file_sim, "\t%.3f", RMSD);
	  }
	}
	
	// Print divergence
	int homo; double qinf=0;
	double DS=Divergence(SI_all[a][i][j], S0, DMAX),
	  DT=Divergence(TM_all[a][i][j], TM0, DMAX),
	  DC=Compute_Dcont(&qinf,CO_all[a][i][j],
			   pi->len,pj->len, &homo,NORM);
	double PC0=(PC0_1+PC_load[3]*qinf)/PC_norm_div;
	double PC_a=PC_all[a][i][j];
	if(DIV_ALIGNED==0){PC_a=(PC_a*PC_norm-na_all[a][i][j])/PC_norm_div;}
	double DP=Divergence(PC_a, PC0, DMAX);
	/*if(file_div_all){
	  fprintf(file_div_all, "\t%.3f\t%.3f\t%.3f\t%.3f", DS, DC, DT, DP);
	  }*/
	if(it){
	  Seq_Id[i][j]=SI_all[a][i][j];
	  Cont_Div[i][j]=DC;
	  TN_Div[i][j]=DS;
	  TM_Div[i][j]=DT;
	  PC_Div[i][j]=DP;
	  // Symmetrize
	  Seq_Id[j][i]=Seq_Id[i][j];
	  Cont_Div[j][i]=DC;
	  TN_Div[j][i]=DS;
	  TM_Div[j][i]=DT;
	  PC_Div[j][i]=DP;
	}		
      } // end loop over alignments
      if(file_sim){
	if(fun_sim_Seq)fprintf(file_sim,"\t%.3f",fun_sim_Seq[i][j]);
	fprintf(file_sim,"\n");
      }
      //if(file_div_all)fprintf(file_div_all,"\n");
      
    } // end j
  } // end i
  
  if(file_sim){
    printf("Similarities between sequences written in %s\n",name_sim);
    fclose(file_sim);
  }
  /*if(file_div_all){
    printf("Divergences between sequences written in %s\n",name_div_all);
    fclose(file_div_all);
    }*/
    
  /*****************************************************************
                 Print statistics for all types
  ******************************************************************/
  if(ALI_SS)printf("%s", ss_def);
  printf("Optimal alignment: %s\n", ali_name[it_opt]);
  
  if(PRINT_ID ){
    char name_id[200];
    sprintf(name_id, "%s.id", name_PC);
    printf("Probabilities of identity written in %s for all pairs\n",name_id);
    FILE *file_id=fopen(name_id, "w");
    fprintf(file_id,"# Conservation properties of aligned residues\n");
    Summary_identical(file_id, it_opt);
    fclose(file_id);
    
    for(i=0; i<ATYPE; i++){
      free(id_aa[i]);
      free(id_sup[i]);
      free(neigh_ali[i]);
      free(neigh_noali[i]);
      free(neigh_noali_aaid[i]);
      free(shift[i]);
      free(id_cont_ct[i]);
      free(ali_cont_ct[i]);
      free(id_cont_sup[i]);
      free(ali_cont_sup[i]);
      free(id_cont_aaid[i]);
      free(ali_cont_aaid[i]);
    }
  }
  
  /******************* Clean *********************************/
  
  for(i=0; i<N_seq; i++){
    int n; if(ALL_PAIRS){n=N_seq;}else{n=i; if(n==0)continue;}
    //printf("Empting Ali_pair[%d] %xd\n", i, Ali_pair[i]);
    //if(Ali_pair[i]){Empty_matrix_i(Ali_pair[i], n);}
  }
  //Empty_matrix_f(PC_Div_opt, N_seq);
  Empty_matrix_i(msa, N_seq);
  Empty_matrix_i(msa_opt, N_seq);
  //Empty_matrix_f(Seq_diff, N_seq);
  if(d2)Empty_matrix_f(d2, L_msa);
  
  /************* Exit if clock violations are not needed *****************/
  if(PRINT_CV==0){
    //printf("Clock violation not required, exiting the program\n");
    printf("Computations finished\n");
    return(0);
  }
  
  /*********************** Clock violations ***********************/
  int Nd=4, dist; // types of distances
  char *name_dist[Nd]; float **Div[Nd], **Div_PDB[Nd];
  for(dist=0; dist<Nd; dist++)name_dist[dist]=malloc(80*sizeof(char));
  name_dist[0]="Tajima-Nei"; Div[0]=TN_Div; 
  name_dist[1]="Cont_Div";   Div[1]=Cont_Div;
  name_dist[2]="TM_Div";     Div[2]=TM_Div;   
  name_dist[3]="PC_Div";     Div[2]=PC_Div_opt;   
  
  for(int it=0; it<RTYPE; it++){
    if(na_all[it])Empty_matrix_f(na_all[it], N_seq);
    if(SI_all[it])Empty_matrix_f(SI_all[it], N_seq);
    if(TM_all[it])Empty_matrix_f(TM_all[it], N_seq);
    if(CO_all[it])Empty_matrix_f(CO_all[it], N_seq);
    if(PC_all[it])Empty_matrix_f(PC_all[it], N_seq);
  }
  
  
  // Neighbor Joining and outgroups
  /* struct treenode *nodes=Neighbor_Joining(TajNei_Div, N_seq);
     int **N_out, ***outgroup=Build_outgroups(&N_out, nodes, N_seq); */
  float **Out_Div=PC_Div_opt;
  if(strcmp(OUTG,"PC")==0){Out_Div=PC_Div_opt;}
  else if(strcmp(OUTG,"CD")==0){Out_Div=Cont_Div;}
  else if(strcmp(OUTG,"TM")==0){Out_Div=TM_Div;}
  else if(strcmp(OUTG,"TN")==0){Out_Div=TN_Div;}
  
  printf("Determining outgroups through neighbor joining ");
  printf("with divergence %s\n", OUTG);
  int **N_out, ***outgroup=Outgroups_NJ(&N_out, Out_Div, N_seq);
  
  // Prepare output
  char head2[5000], ALI[100];
  sprintf(head2,
	  "# %d input proteins %d PDB files\n"
	  "# %d groups of sequences with > %d intragroup seq.id\n"
	  "# Average number of structures per group: %.1f\n"
	  "# Outgroups assigned with NJ using divergence %s "
	  " constructed with alignment %s\n"
	  "# Multiple sequence alignment: %s\n", 
	  N_prot, N_pdb, N_seq, seq_thr, N_pdb/(float)N_seq, file_ali,
	  OUTG, ali_name[o]);
  strcpy(ALI,"SqAli");
  
  char name_out[100];
  Change_ext(name_out, name_in, EXT_CV);
  FILE *file_out=fopen(name_out, "w");
  fprintf(file_out, "# Output of the program %s\n", argv[0]);
  fprintf(file_out, "%s", head);
  fprintf(file_out, "# Prot1, Prot2, number_of_outgroups_C\n");
  fprintf(file_out, "# For each divergence measure d, plot:\n");
  fprintf(file_out, "# d(A,B)\n");
  fprintf(file_out, "# CV(A,B)=sum_C(d(A,C)-d(B,C))/nC*d(A,B)\n"); //^%.2f,alpha
  fprintf(file_out, "# t=|CV(A,B)|/S.E.M.(CV) (Standard Error of Mean)\n");
  fprintf(file_out, "# n_TI number of violations of Triangle Inequality\n");
  fprintf(file_out, "# number of outgroups with minority sign\n");
  fprintf(file_out, "# Note: used outgroups = num_out-num_TI\n");
  fprintf(file_out, "#Prot1 Prot2 L1-L2 num_out");
  int k=5;
  for(dist=0; dist<Nd; dist++){
    fprintf(file_out, " %s:", name_dist[dist]);
    fprintf(file_out, " %d=d %d=CV %d=t_CV %d=n_TI %d=n_sign",
	    k, k+1, k+2, k+3, k+4); k+=5;
  }
  fprintf(file_out, "\n");
  
  // Compute CV
  float **CV_dist[Nd], **t_dist[Nd];
  int **nout_dist[Nd], **nsign_dist[Nd], **nTIV_dist[Nd];
  if(PRINT_AVE){
    for(dist=0; dist<Nd; dist++){
      CV_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      t_dist[dist]=Allocate_mat2_f(N_seq, N_seq);
      nout_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nTIV_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
      nsign_dist[dist]=Allocate_mat2_i(N_seq, N_seq);
    }
  }
  
  int kgap_max=40, bingap=5;
  int **N_gaps=NULL, **All_gaps=NULL, **TIV_gaps=NULL;
  if(PRINT_GAP){
    N_gaps=Allocate_mat2_i(N_seq, N_seq);
    All_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    TIV_gaps=Allocate_mat2_i(Nd, kgap_max+1);
    for(i=0; i<N_seq; i++){
      int i1=i_inp[rep_str[i]];
      for(j=0; j<i; j++){
	int j1=i_inp[rep_str[j]];
	N_gaps[i][j]=Count_gaps(Prot_in[i1].seq, Prot_in[j1].seq, N_ali);
	N_gaps[j][i]=N_gaps[i][j];
      }
    }
  }
  
  int SI_low=0, No_out[Nd]; long Num_out=0;
  for(i=0; i<Nd; i++)No_out[i]=0;
  
  // Sum over pairs of sequences
  for(i=0; i<N_seq; i++){
    for(j=0; j<i; j++){
      Num_out+=N_out[i][j];
      if(Seq_Id[i][j]>SI_thr){SI_low++; continue;}
      fprintf(file_out, "%s\t%s\t%d\t%d",
	      prots[rep_str[i]].domname, prots[rep_str[j]].domname,
	      L_seq[rep_str[i]]-L_seq[rep_str[j]], 
	      N_out[i][j]);
      
      int *outg=outgroup[i][j], used[N_seq];
      for(dist=0; dist<Nd; dist++){
	float **diver=Div[dist], d=diver[i][j];
	// outgroups, triangle inequality, diff.sign
	int nout=0, nTIV=0, nplus=0, kgap=-1, TIV=0; 
	double CV1=0, CV2=0, t=0, nindep=0; // norm=pow(d,alpha);
	for(int k1=0; k1<N_out[i][j]; k1++){
	  // Eliminate outgroups that are very far away
	  int k=outg[k1]; float CV;
	  if((Seq_Id[i][k]<SI_thr)||(Seq_Id[j][k]<SI_thr))continue;
	  CV=diver[i][k]-diver[j][k];
	  if((CV>d)||(CV<-d)){ // check triangle inequality TI
	    used[k1]=0; TIV=1; nTIV++; 
	  }else{
	    used[k1]=1; if(TIV)TIV=0;
	    CV1+=CV; CV2+=CV*CV; nout++; if(CV>0)nplus++;
	    float s_max=0, *Sk=Seq_Id[k];
	    for(int k2=0; k2<k1; k2++){
	      if((used[k2])&&(Sk[outg[k2]]>s_max))s_max=Sk[outg[k2]];
	    }
	    nindep+=s_max; // Sequence identity
	  }
	  if(N_gaps){
	    kgap=N_gaps[i][j]+N_gaps[i][k]+N_gaps[j][k];
	    kgap/=bingap; if(kgap>kgap_max)kgap=kgap_max;
	    All_gaps[dist][kgap]++;
	    if(TIV)TIV_gaps[dist][kgap]++;
	  }
	} // end sum over outgroups
	
	if(nout){CV1/=nout;}
	else{No_out[dist]++;}
	nindep=nout-nindep;
	if(nout <= 1){t=1;}
	else{
	  CV2=(CV2-nout*CV1*CV1)/(nout-1);
	  if(CV2<=0){t=1;}
	  else{t=fabs(CV1)/sqrt(CV2/nindep);}
	}
	CV1/=d;
	
	int nm=nout-nplus; if(nm<nplus)nplus=nm;
	fprintf(file_out, "\t%.3f\t%.3f\t%.1f\t%d\t%d",
		d, CV1, t, nTIV, nplus); //norm
	if(PRINT_AVE){
	  CV_dist[dist][i][j]=CV1;
	  t_dist[dist][i][j]=t;
	  nout_dist[dist][i][j]=nout;
	  nsign_dist[dist][i][j]=nplus;
	  nTIV_dist[dist][i][j]=nTIV;
	}
      } // end dists
      fprintf(file_out, "\n");
    }
  } // end pairs
  fclose(file_out);
  
  if(PRINT_AVE){
    char tmp[200];
    sprintf(tmp, 
	    "# Total number of outgroups: %ld\n# %d pairs with SI < %.2f\n",
	    Num_out, SI_low, SI_thr);
    strcat(head, tmp);
    if(Nd){
      strcat(head, "# pairs with zero outgroups: ");
      for(i=0; i<Nd; i++){
	sprintf(tmp, " %d (%s)", No_out[i], name_dist[i]);
	strcat(head, tmp);
      }
      strcat(head, "%s\n");
    }
    
    for(dist=0; dist<Nd; dist++){
      CV_statistics(name_in, head, OUTG, ALI,
		    name_dist[dist], Div[dist], CV_dist[dist],
		    t_dist[dist], nout_dist[dist],
		    nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		    NULL, 0.00, 1.00, N_seq, dist, "AllFun");
      if(fun_sim_Seq){
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		      fun_sim_Seq, fun_high, 1.00, N_seq, dist, "SameFun");
	CV_statistics(name_in, head, OUTG, ALI,
		      name_dist[dist], Div[dist], CV_dist[dist],
		      t_dist[dist], nout_dist[dist],
		      nsign_dist[dist], nTIV_dist[dist], Seq_Id, SI_thr,
		      fun_sim_Seq, 0.00, fun_low, N_seq, dist, "DiffFun");
      }
      
      Empty_matrix_f(CV_dist[dist], N_seq);
      Empty_matrix_f(t_dist[dist], N_seq);
      Empty_matrix_i(nout_dist[dist], N_seq);
      Empty_matrix_i(nTIV_dist[dist], N_seq);
      Empty_matrix_i(nsign_dist[dist], N_seq);
    }
  }
  if(fun_sim_Seq)Empty_matrix_f(fun_sim_Seq, N_seq);
  
  /***************************************************************************/
  /* Relation between gaps and triangle inequality */
  if(PRINT_GAP){
    Change_ext(name_out, name_in, ".gaps");
    file_out=fopen(name_out, "w");
    for(int dist=0; dist<Nd; dist++){
      fprintf(file_out, "# dist=%s\n", name_dist[dist]);
      fprintf(file_out, "#ngap P(TIV) s.e. num\n");
      double norm=0;
      for(i=0; i<=kgap_max; i++)norm+=All_gaps[dist][i];
      for(i=0; i<=kgap_max; i++){
	if(All_gaps[dist][i]==0)continue;
	float p=(float)TIV_gaps[dist][i]/All_gaps[dist][i];
	fprintf(file_out, "%.1f %.3f %.3f %.3f\n", (i+0.5)*bingap,
		p, sqrt(p*(1-p)/All_gaps[dist][i]),
		All_gaps[dist][i]/norm);
      }
    }
    fclose(file_out);
    printf("Writing %s\n", name_out);
  }
  
  /***************************************************************************/
  return(0);
  }
  //}

void help(char *pname){
  printf("\n\nhelp of program %s\n", pname);
  printf("Author Ugo Bastolla Centro de Biologia Molecular Severo Ochoa "
	 "(CSIC-UAM), Madrid, Spain\n"
	 "Email: <ubastolla@cbm.csic.es>\n\n"
	 "PC_ali performs hybrid multiple structure and sequence alignmentsbased on the structure+sequence similarity score PC_sim, prints pairwise similarity scores and divergence scores and neighbor-joining phylogenetic tree obtained with the hybrid evolutionary divergence measure based on PC_sim. Optionally, it computes violations of the molecular clock for each pair of proteins.\n\n"
	 "It takes as input either a list of PDB files (1=PDB code, 2=chain index, 3=domain decomposition as 1-103/104-207 with SEQRES indexes, 4=folder 5=domain name) in the folder specified by -pdbdir and, if necessary, -pdbdir2, or not aligned sequences (option -seq) or MSA (option -ali). In the last cases PDB codes and names must be specified as sequence names\n\n"
	 "It includes a modification of the needlemanwunsch aligner programmed by Dr. Andrew C. R. Martin in the Profit suite of programs, (c) SciTech Software 1993-2007\n\n"
	 "USAGE:\n"
	 "PC_ali\n"
	 "\t -pdblist <list of PDB files>"
	 "  Format: file_name chain domain dom_name dir"
	 " (only necessary file_name)\n"
	 "\t -pbdir <directory of pdb files>  (default: current directory)\n"
	 "\t -pbdir2 <2nd directory of pdb files>  (default: none)\n"
	 "\t -pdbext <extension of pdb files>  (default: .pdb)\n"
	 "\t -chain_num ! Read chains as numbers instead of char\n\n"
	 "\t -ali <MSA file in FASTA format, with names of PDB files>\n"
	 "\t -seq <sequences in FASTA format, with names of PDB files>\n"
	 "\t # The pdb code is optionally followed by chain index or number\n"
	 "\t # Ex: >1opd.pdb A or >1opd 1 or >1opdA >1opd_A\n"
	 "\t -seq_thr <threshold above which sequences are joined>\n"
	 "\t # For avoiding to join, set -seq_thr 1\n"
	 "\t -id Print propensities between conservation measures\n"
	 "\t -clique ! Initial multiple alignment is based on cliques\n\n");
  printf("Computed similarity measures:\n"
	 "(1) Aligned fraction ali,\n"
	 "(2) Sequence identity SI,\n"
	 "(3) Contact overlap CO,\n"
	 "(4) TM-score TM (Zhang & Skolnick Proteins 2004 57:702)\n"
	 "(5) PC_sim, based on the main Principal Component of the four above similarity scores\n"
	 "They are printed in <>.prot.sim for all pairs of protein sequences, and also for multiple conformations of the same sequence (if present) if required with -print_sim\n\n");

  printf("Computed divergence measures:\n"
	 "(1) Tajima-Nei divergence TN=-log((SI-S0)/(1-S0)) with S0=0.06 (Tajima F & Nei 1984, Mol Biol Evol 1:269),\n"
	 "(2) Contact_divergence CD=-log((q-q0(L))/(1-q0(L)) (Pascual-Garcia et al Proteins 2010 78:181-96),\n"
	 "(3) TM_divergence=-log((TM-TM0)/(1-TM0)), TM0=0.167.\n"
	 "(4) PC_divergence=-log((PC-PC0)/(1-PC0)), PC0 linear combination of S0, TM0, CO(L) and nali0=0.5.\n"
	 "They are printed in <>.prot.div for all pairs of protein sequences, and also for multiple conformations of the same sequence (if present) if required with -print_div\n\n");

  printf("Flux of the program:\n"
	 "(1) In the modality -ali, the program starts from the pairwise alignments obtained from the input MSA. In the modality -seq the starting pairwise alignments are built internally.\n"
	 "(2) The program then modifies the pairwise alignments by targeting PC_sim. The similarity matrix is constructed recursively, using the input pairwise alignment for computing the shared contacts and the distance after optimal superimposition (maximizing the TM score) for all pairs of residues and obtaining a new alignments. Two iterations are usually enough for getting good results. Optionally, for the sake of comparison, the program can target the TM score (-ali_tm), the Contact Overlap (-ali_co) and the secondary structure superposition (-ali_ss).\n"
	 "(3) Sequences with > seq_thr identity are joined together, to accelerate computations and reduce the output size. They represent different conformations of the same protein. The structural similarity (divergence) between two proteins is computed as the maximum (minimum) across all the examined conformations. The clustering can be avoided by setting -seq_thr 1\n"
	 "(4) Then, the program builds progressive multiple alignments. If the option -clique is set, the starting MSA is based on the maximal cliques of the pairwise alignments, which does not require neither a guide tree nor gap penalty parameters, however it can become slow for large data sets.\n"
	 "(5) Finally, the program runs iteratively progressive multiple alignments using as guide tree the average linkage tree obtained with the PC_Div divergence measure of the previous step and using as starting alignment the previous multiple alignment. The best MSA is selected as the one with the maximum value of the average PC similarity score.\n"
	 "(6) The program prints the optimal MSA and the Neighbor Joining tree obtained from the corresponding PC_Div divergence measure.\n"
	 "(7) If the options -print_sim or -print_div are set, the program prints in files <>.prot.sim and <>.prot.div similarity and divergence scores for the input MSA (if present) and for the final MSA.\n"
	 "(8) If -print_pdb is set, the program prints the multiple superimposition obtained by maximizing the TM score\n"
	 "(9) Furthermore, if -print_cv is set, the program computes and prints for all four divergence measures the violations of the molecular clock averaged over all possible outgroups identified with the Neighbor-Joining criterion, and the corresponding significance score.\n\n");

  printf("COMPILE:\n"
	 ">unzip PC_ali.zip\n"
	 ">make\n"
	 ">cp PC_ali ~/bin/ (or whatever path directory you like)\n"
	 "\n"
	 "RUN:\n"
	 "1) with pdblist:\n"
	 " >PC_ali -pdblist 1.10.287.110.SI60..pdblist"
	 " -pdbdir <PDBPATH>\n"
	 "(all files listed in 1.10.287.110.SI60.pdblist must be in "
	 "current folder or in pdbdir)\n\n"
	 "2) with not aligned sequenses:\n"
	 ">PC_ali -seq 50044_Mammoth.aln -pdbdir <PDBPATH>\n"
	 "(all PDB files named in 50044_Mammoth.aln must be in <PDBPATH>)\n\n"
	 "3) with aligned sequenses:\n"
	 ">PC_ali -ali 50044_Mammoth.aln -pdbdir <PDBPATH>\n"
	 "(all PDB files named in 50044_Mammoth.aln must be in <PDBPATH>)\n\n"
	 );

  printf("OUTPUT (for each pair of proteins):\n"
	 "-------\n"
	 "MSA (.msa),\n"
	 "NJ tree (.tree), \n"
	 "structure similarity (.sim) and structure divergence scores (.div) for each protein pair,\n"
	 "correlations between different types of sequence and structure identity (.id), MSA of secondary structure (.ss.msa)\n"
	 "-------\n");
  
  printf("==========================================================\n\n"
	 "Inputs and Options:\n"
	 "\t -pdblist <list of PDB files + chain labels + dir (1 or 2)>\n"
	 "\t -seq <sequences in FASTA format, with names of PDB files>\n"
	 "\t -ali <MSA file in FASTA format, with names of PDB files>\n"
	 "\t# The pdb code is optionally followed by the chain index\n"
	 "\t# Ex: >1opd.pdb A or >1opdA or >1opd_A\n\n"
	 "\t -pdbdir <directory of pdb files>  (default: current directory)\n"
	 "\t -pdbdir2 <2nd directory of pdb files>  (default: current)\n"
	 "\t -pdbext <extension of pdb files>  (default: .pdb)\n"
	 "#### Optional parameters:\n"
	 "\t -out <Name of output files> (default: alignment file)\n"
	 "\t -seq_thr     ! Identity above which sequences are joined "
	 "(def. %.3f)\n"
	 "\t -print_pdb    ! Print multiple structure superimposition\n"
	 "\t -print_sim    ! Print similarity measures for all pairs\n"
	 "\t -print_div    ! Print PC divergence for all pairs\n"
	 "\t -id           ! Print propensities between conservation measures\n"
	 "\t -clique ! Initial alignment is based on cliques\n"
	 "\t -ali_tm     ! Make pairwise alignments that target TM score\n"
	 "\t -ali_co     ! Make pairwise alignments that target Cont Overlap\n"
	 "\t -ali_ss     ! Make alignments that target sec.structure\n"
	 "\t -ss_mult    ! target sec.structure with multiple alignment\n"
	 "\t -shift_max <Maximum shift for targeting sec.str.>\n"
	 "\t -id     ! Print statistics of identical residues\n"
	 "\t -print_cv     ! Print clock violations\n"
	 "\t -print_fasta  ! Print sequences of selected domains\n"
	 "\t -exit_fasta  ! Exit program after printing sequences\n"
	 "\t -func <file with function similarity for pairs of proteins>\n"
	 "\n", seq_thr);
  exit(8);
}

int Remove_gap_cols(int **msa, int N_seq, int L_msa){
  int L_ali=0, jcol[L_msa], i, j;
  for(j=0; j<L_msa; j++){
    int keep=0;
    for(i=0; i<N_seq; i++){if(msa[i][j]>=0){keep=1; break;}}
    if(keep){jcol[j]=L_ali; L_ali++;}
    else{jcol[j]=-1;}
  }
  if(L_ali!=L_msa){
    for(j=0; j<L_msa; j++){
      int k=jcol[j]; if(k<0 || k==j){continue;}
      // Copy column j at place k
      for(i=0; i<N_seq; i++){msa[i][k]=msa[i][j];}
    }
  }
  return(L_ali);
}

float Mean_freq(double *se, double *relerr2, double sum, double tot){
  if(tot==0)return(0);
  float p=sum/tot;
  if(sum){
    *relerr2=((1-p)/sum); //(sum*sum)
    *se=p*sqrt(*relerr2);
  }else{
    *relerr2=1; *se=0;
  }
  return(p);
}

void Print_ali_ss(int *ali_i, short *ss3_i, int *ali_j, short *ss3_j,
		  int N_ali)
{
  for(int i=0; i<N_ali; i++){
    if(ali_i[i]<0){if(ali_j[i]>=0)printf("-");}
    else{printf("%c", SS_code[ss3_i[ali_i[i]]]);}
  }
  printf("\n");
  for(int i=0; i<N_ali; i++){
    if(ali_j[i]<0){if(ali_i[i]>=0)printf("-");}
    else{printf("%c", SS_code[ss3_j[ali_j[i]]]);}
  }
  printf("\n");
}

int Get_pdb_list(struct Prot_input **Prot_input, char *input)
{
  // Format of PDB list: PDB chain domain dir or PDB chain dir domain?
  // What about domain name? What if one of them is missing?
  // Format a 5 col : PDB chain domain domname dir. Missing columns should be NA
  // Format b 3 col:  PDB chain dir
  // Example from Georgia: 1fdr_	A	2	0.999

  // Domain example: 10mhA01	1-186,285-327
  // Open file
  if(input[0]=='\0')return(0);
  FILE *file_in=fopen(input, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", input); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]!='#')n++;
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no proteins found in file %s\n", input); exit(8);
  }
  printf("%d proteins found in %s\n", n, input);

  // Allocate and read PDB codes
  file_in=fopen(input, "r");
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  n=0; int isdom=-2;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    int dir=1;
    char name[80], chain[10], word3[200], domname[40]; 
    int arg=sscanf(string, "%s%s%s%s%d", name, chain, word3, domname, &dir);
    struct Prot_input *Pn=(*Prot_input)+n;
    strcpy(Pn->code, name);
    Pn->chain=' ';
    Pn->chain_num=-1;
    Pn->domname[0]='\0';
    Pn->domain[0]='\0';
    Pn->nfrag=1;
    Pn->ini_frag=NULL;
    Pn->end_frag=NULL;

    // get chain
    if(CHAIN_NUM==0){
      if((arg>1)&&(chain[0]!='\0')&&(chain[0]!='\n')){
	Pn->chain=chain[0];
      }
    }else{
      if(arg>1){sscanf(chain, "%d", &(Pn->chain_num));}
    }

    // Read domain
    if(arg>=4){strcpy(Pn->domname, domname);}
    else{sprintf(Pn->domname, "%s%c",Pn->code, Pn->chain);}
    printf("%s", Pn->domname);

    printf(" pdb=%s chain=%c %d", Pn->code, Pn->chain, Pn->chain_num);

    // What is third column?
    if(arg>=3 && isdom==-2){
      // if contains - it is domain, if 1 or 2 it is dir
      char *ptr=word3;
      while(*ptr!='\0'){if(*ptr=='-'){isdom=1; break;} ptr++;}
      if(isdom<0){
	sscanf(word3, "%d", &dir); if(dir==1 || dir==2){isdom=0;}
      }
    }

    // Read domain or dir
    if(arg>=3){
      if(isdom==1){
	strcpy(Pn->domain, word3);
	Pn->nfrag=Read_domain(&(Pn->ini_frag), &(Pn->end_frag), Pn->domain);
	printf(" %d domains: %s\n", Pn->nfrag, Pn->domain);
      }else if(isdom==0){
	sscanf(word3, "%d", &dir);
	printf(" folder=%d\n", dir);
      }else{
	printf("\nERROR, col 3 must be domain (ex. 1-78) or folder (1 or 2)\n"
	       "isdom=%d\nJust read in %s:\n %s", isdom, input, string);
	exit(8);
      }
    }

    // Check dir
    if(dir==2 && PDBDIR2[0]=='\0'){
      printf("ERROR reading %s, dir (col 3 or 5) is 2 but "
	     "PDBDIR2 not given\n",
	     "Just read: %s", input, string); exit(8);
    }else if(dir!=1 && dir!=2){
      printf("ERROR reading %s, dir (col 3 or 5) must be either 1 or 2\n"
	     "Just read: %s", input, string); exit(8);
    }
    Pn->dir=dir-1;

    (*Prot_input)[n].seq=NULL;
    (*Prot_input)[n].len=0;
    n++;
  }
  fclose(file_in);

  return(n);
}

int Read_domain(int **ini_frag, int **end_frag, char *domain)
{
  char dom_store[200]; strcpy(dom_store, domain);
  char *frag=dom_store; int n_frag=1, m_frag=0, isfrag=0;
  while(*frag!='\0'){
    if(*frag=='-'){m_frag++; isfrag=1;}
    else if(*frag==','){n_frag++;}
    frag++;
  }
  if(isfrag==0){
    printf("ERROR, %s wrong domain format (no -)\n", domain); exit(8);
  }
  if(n_frag!=m_frag){
    printf("ERROR, %s wrong domain format (%d - but %d ,)\n",
	   domain, m_frag, n_frag); exit(8);
  }
  *ini_frag=malloc(n_frag*sizeof(int));
  *end_frag=malloc(n_frag*sizeof(int));
  frag=dom_store; char *frag1=frag; n_frag=0;
  while(*frag !='\0'){
    if(*frag=='-'){
      *frag=' '; sscanf(frag1, "%d", (*ini_frag)+n_frag); *frag='-';
      frag1=frag+1;
    }else if(*frag==','){
      *frag=' '; sscanf(frag1, "%d", (*end_frag)+n_frag);  *frag='-';
      frag1=frag+1;
      n_frag++;
    }
    frag++;
  }
  sscanf(frag1, "%d", (*end_frag)+n_frag); n_frag++;
  return(n_frag);
}

int Get_sequences(struct Prot_input **Prot_input, int *Nali, char *file_ali,
		  int INP_MSA) //char *PDBDIR, char *PDB_EXT, 
{
  // Open file
  if(file_ali[0]=='\0')return(0);
  FILE *file_in=fopen(file_ali, "r");
  if(file_in==NULL){
    printf("ERROR, alignment file %s does not exist\n", file_ali); exit(8);
  }
  // Count proteins and read path
  char string[1000]; int n=0;
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='>')n++;
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no sequence found in file %s\n", file_ali); exit(8);
  }
  printf("%d sequences found in %s\n", n, file_ali);

  int lali[n];

  // Allocate and read sequences
  int l=0, i;
  //int LMAX=10000; char Seq[LMAX];
  *Prot_input=malloc(n*sizeof(struct Prot_input));
  struct Prot_input *Pi=*Prot_input;
  for(i=0; i<n; i++){
    Pi->nfrag=1;
    Pi->ini_frag=NULL;
    Pi->end_frag=NULL;
    Pi->domname[0]='\0';
    Pi->domain[0]='\0';
    Pi->dir=0;
    Pi++;
  }
  n=-1; *Nali=0;

  char word1[40], word2[40], word3[40], *s=NULL; int d1;
  file_in=fopen(file_ali, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    /*if(strncmp(string,"PDBDIR",6)==0){
      if(PDBDIR[0]!='\0'){
	printf("WARNING reading line %sPDBDIR is already set to %s "
	       "keeping it\n", string, PDBDIR);
      }else{
	sscanf(string+7, "%s", PDBDIR);
      }
      continue;
      }*/
    if(string[0]=='>'){
      n++;
      (*Prot_input)[n].chain=' ';
      (*Prot_input)[n].dir=0;
      strcpy(word1, "");
      strcpy(word2, "");
      int c=sscanf(string+1, "%s%s%s%d",
		   (*Prot_input)[n].domname, word1, word2, &d1);
      char *name=(*Prot_input)[n].domname;
      int nc=0; char *ch=word1; while(*ch!='\0'){ch++; nc++;}

      // Get chain
      if(CHAIN_NUM==0){
	if(c==1){
	  if(name[4]=='_'){(*Prot_input)[n].chain=name[5];}
	  else{(*Prot_input)[n].chain=name[4];}
	}else if(nc==1){
	  (*Prot_input)[n].chain=word1[0];
	}
      }else{
	if(c==1){
	  if(name[4]=='_'){
	    sscanf(name+5, "%d", &((*Prot_input)[n].chain_num));
	  }
	}else{
	  sscanf(word1, "%d", &((*Prot_input)[n].chain_num));
	}
      }

      // Get domain
      if(c>2){
	strcpy((*Prot_input)[n].domain, word2);
      }

      // Get code
      name[4]=' ';
      sscanf(name, "%s", (*Prot_input)[n].code);

      // Get dir
      if(c>3){(*Prot_input)[n].dir=d1;}

      printf("%s %c\n", (*Prot_input)[n].code, (*Prot_input)[n].chain);
      if(INP_MSA && (*Nali==0)&&(l)){*Nali=l;} l=0;

      /*if((*Nali==0)&&(l)){
	*Nali=l;
	(*Prot_input)[0].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[0].seq;
	for(l=0; l<*Nali; l++){*s=Seq[l]; s++;}
      }
      if(*Nali){
	(*Prot_input)[n].seq=malloc(*Nali*sizeof(char));
	s=(*Prot_input)[n].seq; l=0;
      }else{
	s=Seq; l=0;
	} */
    }else{
      char *c=string;
      while((*c!='\n')&&(*c!='\r')&&(*c!='\0')&&(*c!=' ')){
	/*if(l > LMAX){
	  printf("ERROR, alignment %s length %d > maximum allowed %d\n",
		 (*Prot_input)[n].domname, l, LMAX);
	  printf("Increase LMAX in code %s\n", CODENAME); exit(8);
	  }*/
	if(INP_MSA && (l>*Nali)&&(*Nali)){
	  // Check if all additional characters are gaps
	  int allgaps=1; char *d=c; int i=l;
	  while((*d!='\n')&&(*d!='\r')&&(*d!='\0')&&(*d!=' ')){
	    if(*d!='-'){allgaps=0;} d++; i++;
	  }
	  if(allgaps==0){
	    printf("ERROR, too many columns in sequence %s",
		   (*Prot_input)[n].domname);
	    printf(" Expected %d, found %d\n", *Nali, i);
	    for(i=0; i<l; i++)printf("%c", (*Prot_input)[n].seq[i]);
	    printf("+");
	    char *d=c; int k=n-1;
	    while((*d!='\n')&&(*d!='\r')&&(*d!='\0')&&(*d!=' ')){
	      printf("%c",*d); d++;
	    }
	    printf("\n");
	    printf("Previous seq %s:\n",(*Prot_input)[k].domname);
	    for(i=0; i<*Nali; i++)printf("%c", (*Prot_input)[k].seq[i]);
	    printf("\n");
	    exit(8); // end error
	  }else{ // allgaps==1
	    printf("WARNING, length %d >= Nali=%d but all gaps\n",l,*Nali);
	    break;
	  } // end test all gaps
	}// end l>Nali
	//*s=*c; s++;
	l++; c++;
      }// end read line
      lali[n]=l;
    } // end type of line
  } // end read file
  n++;
  fclose(file_in);
  printf("%d sequences with %d columns found in MSA %s\n",
	 n, *Nali, file_ali);

  for(int k=0; k<n; k++){
      if(lali[k]>(*Nali)){(*Nali)=lali[k];}
  }
  for(int k=0; k<n; k++){
    (*Prot_input)[k].seq=malloc(*Nali*sizeof(char));
    if(lali[k]<(*Nali)){ // add gaps
      char *s=(*Prot_input)[k].seq+lali[k];
      for(int i=lali[k]; i<(*Nali); i++){*s='-'; s++;}
    }
  }
  n=-1;

  file_in=fopen(file_ali, "r");
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#' || strncmp(string,"PDBDIR",6)==0)continue;
    if(string[0]=='>'){
      n++;  s=(*Prot_input)[n].seq; l=0;
    }else{
      char *c=string;
      while((*c!='\n')&&(*c!='\r')&&(*c!='\0')&&(*c!=' ')){
	*s=*c; s++; l++; c++;
      }// end read line
    }
  }
  n++;
  fclose(file_in);

  /*if(PDBDIR[0]=='\0'){
    printf("WARNING , PDBDIR not set, using current directory\n");
    strcpy(PDBDIR, "./");
    }*/

  return(n);
}

void Write_identity(int type,
		    struct protein *proti, struct protein *protj,
		    int *ali_ij, int *id_aa, int *id_sup, int *shift,
		    int *neigh_ali,int *neigh_noali, int *neigh_noali_aaid,
		    int *ali_cont_ct, int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid, int *id_cont_aaid)
{
  if(ini_sum==0){
    ini_sum=1;
    for(int i=0; i<CTYPE; i++){
      sum_all_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_ali_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_aaid_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_sup_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_aaid_sup_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_noid_nosup_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_ali_cont_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_id_cont_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_neigh_ali_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_neigh_noali_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_neigh_noali_aaid_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_sup_noneigh_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_shift_sup_noneigh_ct[i]=Allocate_mat2_d(2, IBIN1);
      sum_ali_cont_sup[i]=Allocate_mat2_d(3, IBIN1);
      sum_id_cont_sup[i]=Allocate_mat2_d(3, IBIN1);
      sum_ali_cont_aaid[i]=Allocate_mat2_d(3, IBIN1);
      sum_id_cont_aaid[i]=Allocate_mat2_d(3, IBIN1);
      for(int k=0; k<3; k++){
	sum_ali_cont_aaid_sup[i][k]=Allocate_mat2_d(2, IBIN1);
	sum_id_cont_aaid_sup[i][k]=Allocate_mat2_d(2, IBIN1);
      }
    }  
    for(int k=0; k<2; k++){
      sum_SS_SI[k]=0; sum_SS_TM[k]=0; sum_SS_CO[k]=0;
      sum_SS_SI_TM[k]=malloc(2*sizeof(double));
      sum_SS_SI_CO[k]=malloc(2*sizeof(double));
      sum_SS_TM_CO[k]=malloc(2*sizeof(double));
      for(int j=0; j<2; j++){
	sum_SS_SI_TM[k][j]=0; sum_SS_SI_CO[k][j]=0; sum_SS_TM_CO[k][j]=0;
      }
    }
  }

  int k, noali[2], nali[2], n_sup[2], n_aaid[2], n_aaid_sup[2], n_noid_nosup[2],
    n_neigh_ali[2], n_neigh_noali[2], n_neigh_noali_aaid[2],
    n_sup_noneigh[2], shift_sup_noneigh[2];
  for(k=0; k<2; k++){
    noali[k]=0; nali[k]=0; n_aaid[k]=0; n_sup[k]=0;
    n_aaid_sup[k]=0; n_noid_nosup[k]=0; 
    n_neigh_ali[k]=0; n_neigh_noali[k]=0; n_neigh_noali_aaid[k]=0;
    n_sup_noneigh[k]=0; shift_sup_noneigh[k]=0;
  }
  // Not aligned sites
  if(proti->len<protj->len){
    Count_noali(noali, ali_ij,  proti->ncont, proti->len, c_ave);
  }else{
    int ali2[protj->len]; Invert_ali(ali2, protj->len, ali_ij, proti->len);
    Count_noali(noali, ali2, protj->ncont, protj->len, c_ave);
  }
  for(int i=0; i<proti->len; i++){
    if(ali_ij[i]<0){
      continue; // These are not errors in alignments!
      if(neigh_noali[i]==0)continue;
      if(proti->ncont[i]<=c_ave){k=0;}// few contacts
      else{k=1;} // many contacts
      n_neigh_noali[k]++;
      if(neigh_noali_aaid[i])n_neigh_noali_aaid[k]++;
      continue;
    }
    float c=sqrt(proti->ncont[i]*protj->ncont[ali_ij[i]]);
    if(c<=c_ave){k=0;} // few contacts
    else{k=1;}  // many contacts
    nali[k]++;
    if(id_aa[i]){n_aaid[k]++; if(id_sup[i])n_aaid_sup[k]++;}
    else{if(id_sup[i]==0)n_noid_nosup[k]++;}
    if(id_sup[i]){
      n_sup[k]++;
      if(neigh_ali[i]==0){
	n_sup_noneigh[k]++; shift_sup_noneigh[k]+=shift[i];
      }
    }
    if(neigh_ali[i]){n_neigh_ali[k]++;}
    else if(neigh_noali[i]){n_neigh_noali[k]++;
      if(neigh_noali_aaid[i])n_neigh_noali_aaid[k]++;
    }
  }

  //int j=(n_aaid[0]+n_aaid[1])*IBIN/(float)(nali[0]+nali[1]);
  int j=(n_aaid[0]+n_aaid[1]), ntot=nali[0]+nali[1];
  if(ntot){j*=(IBIN/(float)ntot);}
  if(j>IBIN){j=IBIN;}
  for(k=0; k<2; k++){
    sum_all_ct[type][k][j]+=(noali[k]+nali[k]);
    sum_ali_ct[type][k][j]+=nali[k];
    sum_sup_ct[type][k][j]+=n_sup[k];
    sum_aaid_ct[type][k][j]+=n_aaid[k];
    sum_aaid_sup_ct[type][k][j]+=n_aaid_sup[k];
    sum_noid_nosup_ct[type][k][j]+=n_noid_nosup[k];
    sum_neigh_ali_ct[type][k][j]+=n_neigh_ali[k];
    sum_neigh_noali_ct[type][k][j]+=n_neigh_noali[k];
    sum_neigh_noali_aaid_ct[type][k][j]+=n_neigh_noali_aaid[k];
    sum_sup_noneigh_ct[type][k][j]+=n_sup_noneigh[k];
    sum_shift_sup_noneigh_ct[type][k][j]+=shift_sup_noneigh[k];
    //
    sum_ali_cont_ct[type][k][j]+=ali_cont_ct[k];
    sum_id_cont_ct[type][k][j]+=id_cont_ct[k];
  }
  for(k=0; k<3; k++){
    sum_ali_cont_sup[type][k][j]+=ali_cont_sup[k];
    sum_id_cont_sup[type][k][j]+=id_cont_sup[k];
    sum_ali_cont_aaid[type][k][j]+=ali_cont_aaid[k];
    sum_id_cont_aaid[type][k][j]+=id_cont_aaid[k];
  }
}


void Get_file_name(char *name, char *file){
  char *name_nodir=Name_nodir(file);
  Name_noext(name, name_nodir);
}

void Set_contact_type(){

  if(CONT_TYPE=='a'){strcpy(CONT_STRING, "Alpha");}
  else if(CONT_TYPE=='b'){strcpy(CONT_STRING, "Beta");}
  else if(CONT_TYPE=='c'){strcpy(CONT_STRING, "All atoms");}
  else{
    printf("WARNING, undefined contact %c\n", CONT_TYPE);
    CONT_TYPE='c'; strcpy(CONT_STRING, "All atoms");
    printf("Using default %s\n", CONT_STRING);
  }
  // Default type of contacts?
  if((CONT_TYPE!=CONT_TYPE_DEF)||(CONT_THR!=CONT_THR_DEF)||
     (IJ_MIN!=IJ_MIN_DEF))CONT_DEF=0;
}

int Count_AA(char *seq, int N_ali){
  int L=0; char *s=seq;
  for(int i=0; i<N_ali; i++){if(*s!='-')L++; s++;}
  return(L);
}

int Seq_differences(int *id, char *seq1, char *seq2, int N_ali){
  // Consider only aligned positions
  int d=0; (*id)=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){ // ||
      if(*s1!=*s2){d++;}else{(*id)++;}
    }
    s1++; s2++;
  }
  return(d);
}


int Count_gaps(char *seq1, char *seq2, int N_ali){
  int ngap=0, open=0; char *s1=seq1, *s2=seq2;
  for(int i=0; i<N_ali; i++){
    if((*s1!='-')&&(*s2!='-')){open=0;}
    else if(open==0){ngap++; open=1;}
    s1++; s2++;
  }
  return(ngap);
}

float Min_CV(int i, int j, int k, int **conformation, int *N_conf, float **div)
{
  float CV=10000;
  for(int k1=0; k1<N_conf[k]; k1++){
    float *dk=div[conformation[k][k1]];
    for(int i1=0; i1<N_conf[i]; i1++){
      float d_ki=dk[conformation[i][i1]];
      for(int j1=0; j1<N_conf[j]; j1++){
	float cv=d_ki-dk[conformation[j][j1]];
	if(fabs(cv)<fabs(CV)){CV=cv;}
      }
    }
  }
  return(CV);
}

/*
void Change_conformations(int **conformation, int N_seq, int *N_conf,
			  int *ali_str)
{
  for(int i=0; i<N_seq; i++){
    for(int k=0; k<N_conf[i]; k++){
      int l=conformation[i][k];
      conformation[i][k]=ali_str[l];
    }
  }
}
*/
float Min_dist(int i, int j, int **conformation, int *N_conf, float **div)
{
  float d=1000; int ini=1;
  for(int i1=0; i1<N_conf[i]; i1++){
    float *d_i=div[conformation[i][i1]];
    for(int j1=0; j1<N_conf[j]; j1++){
      float d_il=d_i[conformation[j][j1]];
      if(ini){d=d_il; ini=0;}
      else if(d_il<d){d=d_il;}
    }
  }
  return(d);
}

/*
int *Match_alignments(struct Prot_input *Prot1,
		      struct Prot_input *Prot2, int N)
{
  int *ali_str=malloc(N*sizeof(int));
  struct Prot_input *P1=Prot1;
  for(int i=0; i<N; i++){
    ali_str[i]=-1;
    for(int j=0; j<N; j++){
      if((strcmp(Prot2[j].domname, P1->domname)==0)&&
	 (Prot2[j].chain==P1->chain)){
	ali_str[i]=j; break;
      }
    }
    if(ali_str[i]<0){
      printf("WARNING, protein %s%c i=%d N=%d not matched\n",
	     P1->domname, P1->chain, i, N);
      free(ali_str); return(NULL);
    } 
    P1++;
  }

  return(ali_str);
}

*/

char *Name_nodir(char *file){
  char *ini=file, *c=file;
  while(*c !='\0'){if(*c=='/'){ini=c+1;} c++;}
  if(*ini=='\0'){ini=file;}
  return(ini);
}

void Name_noext(char *file_new, char *file){
  strcpy(file_new, file);
  char *c=file_new, *dot=NULL;
  while(*c !='\0'){if(*c=='.'){dot=c;} c++;}
  if(dot!=NULL){c=dot; while(*c !='\0'){*c='\0';}}
  return;
}

void Get_input(char *file_ali, char *file_list, char *file_fun,
	       char *name_in,  char *PDB_DIR, char *PDB_EXT, char *OUTG,
	       int *NORM, int *ALI_SS, int *SHIFT_MAX, int *SS_MULT,
	       int *PRINT_SIM, int *PRINT_CV, int *PRINT_DIV,
	       //int *PRINT_PAIR, int *PRINT_CLIQUE,
	       int argc, char **argv)
{
  printf("Starting %s\n", argv[0]);
  if(argc<2)help(argv[0]);

  char input[80]="";
  for(int i=1; i<argc; i++){
    if(strncmp(argv[i], "-h", 2)==0){help(argv[0]);}
    else if(strncmp(argv[i], "-file", 5)==0){
      i++; strcpy(input, argv[i]); break;
    }
  }

  int ali=0;
  if(input[0]!='\0'){
    FILE *file_in=fopen(argv[1], "r");
    if(file_in==NULL){goto read_command_line;}
    printf("Reading parameters from file %s\n", argv[1]);
    char string[100], read[80];
    while(fgets(string, sizeof(string), file_in)!=NULL){
      if(string[0]=='#')continue;
      if(strncmp(string, "ALI=", 4)==0){
	sscanf(string+4,"%s", file_ali);
	printf("File with multiple sequence alignment: %s\n", file_ali);
	ali=1;
      }else if(strncmp(string, "PDBLIST=", 7)==0){
	sscanf(string+8,"%s", file_list);
	printf("List of PDB files: %s\n", file_list);
      }else if(strncmp(string, "NAME", 4)==0){
	sscanf(string+5, "%s", name_in);
      }else if(strncmp(string, "FUN_SIM", 7)==0){
	sscanf(string+8, "%s", file_fun);
      }else if(strncmp(string, "OUTGROUP", 8)==0){
	sscanf(string+9, "%s", read);
	if((strcmp(read,"TN")!=0) && 
	   (strcmp(read,"CD")!=0) && 
	   (strcmp(read,"TM")!=0)){
	  printf("WARNING, outgroup method %s not implemented, using default\n",
		 read);
	}else{
	  strcpy(OUTG, read);
	}
	printf("Outgroup method set to %s\n", OUTG);
      }else if(strncmp(string, "PDBDIR", 6)==0){
	sscanf(string+7,"%s", PDB_DIR);
	printf("Directory for PDB files: %s\n", PDB_DIR);
      }else if(strncmp(string, "NORM", 4)==0){
	char tmp[80]; sscanf(string+5,"%s", tmp);
	if((strcmp(tmp,"MIN")==0)||
	   (strcmp(tmp,"Min")==0)||
	   (strcmp(tmp,"min")==0)){*NORM=0;}
	else if((strcmp(tmp,"MAX")==0)||
		(strcmp(tmp,"Max")==0)||
		(strcmp(tmp,"max")==0)){*NORM=1;}
	else if((strcmp(tmp,"MEAN")==0)||
		(strcmp(tmp,"Mean")==0)||
		(strcmp(tmp,"mean")==0)){*NORM=2;}
	else{
	  printf("WARNING, %s is not an allowed normalization type\n",tmp);
	  printf("Keeping default MIN\n");
	}
      }else if(strncmp(string, "ALI_TM", 6)==0){
	sscanf(string+7,"%d", &ALI_TM);
      }else if(strncmp(string, "ALI_CO", 6)==0){
	sscanf(string+7,"%d", &ALI_CO);
      }else if(strncmp(string, "ALI_PC", 6)==0){
	sscanf(string+7,"%d", &ALI_PC);
      }else if(strncmp(string, "ALI_SS", 6)==0){
	sscanf(string+7,"%d", ALI_SS);
      }else if(strncmp(string, "SS_MULT", 7)==0){
	sscanf(string+8,"%d", SS_MULT);
      }else if(strncmp(string, "SHIFT_MAX", 9)==0){
	int tmp; sscanf(string+10,"%d", &tmp);
	if(tmp>=0){*SHIFT_MAX=tmp;}
	else{printf("WARNING, %d not allowed value for SHIFT_MAX\n", tmp);}
      }else if(strncmp(string, "PDBEXT", 6)==0){
	sscanf(string+7, "%s", PDB_EXT);
      }else if(strncmp(string, "PRINT_SIM", 9)==0){
	sscanf(string+10, "%d", PRINT_SIM);
      }else if(strncmp(string, "PRINT_CV", 8)==0){
	sscanf(string+9, "%d", PRINT_CV);
      }else if(strncmp(string, "PRINT_DIV", 9)==0){
	sscanf(string+10, "%d", PRINT_DIV);
	/*}else if(strncmp(string, "PRINT_PAIR", 10)==0){
	sscanf(string+11, "%d", PRINT_PAIR);
      }else if(strncmp(string, "PRINT_CLIQUE", 12)==0){
      sscanf(string+13, "%d", PRINT_CLIQUE);*/
      }else{
	printf("WARNING, unknown command %s", string);
      }
    }

    fclose(file_in);
  }

 read_command_line:
  // Read from command line
  for(int i=1; i<argc; i++){
    if(strcmp(argv[i], "-ali")==0){
      i++; sscanf(argv[i], "%s", file_ali);
      if(ali){printf("WARNING, seq already provided\n");}
      ali=1; INP_MSA=1;
    }else if(strcmp(argv[i], "-chain_num")==0){
      CHAIN_NUM=1;
    }else if(strcmp(argv[i], "-pdblist")==0){
      i++; sscanf(argv[i], "%s", file_list);
    }else if(strcmp(argv[i], "-seq")==0){
      i++; 
      if(ali){printf("WARNING, MSA already provided, using it\n");}
      else{sscanf(argv[i], "%s", file_ali); INP_MSA=0;}
      ali=1;
    }else if(strcmp(argv[i], "-out")==0){
      i++; sscanf(argv[i], "%s", name_in);
    }else if(strcmp(argv[i], "-pdbdir2")==0){
      i++; sscanf(argv[i], "%s", PDBDIR2);
      printf("Directory for PDB files: %s\n", PDBDIR2);
    }else if(strcmp(argv[i], "-pdbdir")==0){
      i++; sscanf(argv[i], "%s", PDB_DIR);
      printf("Directory for PDB files: %s\n", PDB_DIR);
    }else if(strcmp(argv[i], "-pdbext")==0){
      i++; sscanf(argv[i], "%s", PDB_EXT);
    }else if(strcmp(argv[i], "-clique")==0){
      MAKE_CLIQUE=1;
      /*}else if(strcmp(argv[i], "-print_clique")==0){
       *PRINT_CLIQUE=1;
       }else if(strcmp(argv[i], "-print_pair")==0){
       *PRINT_PAIR=1;
       }else if(strcmp(argv[i], "-print_msa")==0){
       MAKE_MSA=1;*/
    }else if(strcmp(argv[i], "-seq_thr")==0){
      i++; sscanf(argv[i], "%f", &seq_thr);
      if(seq_thr>1){
	printf("ERROR in input, seq_thr must be <= 1\n");
	exit(8);
      }
    }else if(strcmp(argv[i], "-print_pdb")==0){
      PRINT_PDB=1;
    }else if(strcmp(argv[i], "-print_fasta")==0){
      PRINT_FASTA=1;
    }else if(strcmp(argv[i], "-exit_fasta")==0){
      EXIT_FASTA=1;

    }else if(strcmp(argv[i], "-id")==0){
      PRINT_ID=1;
      printf("Printing propensities between conservation measures\n");
    }else if(strcmp(argv[i], "-print_sim")==0){
      *PRINT_SIM=1;
    }else if(strcmp(argv[i], "-print_cv")==0){
      *PRINT_CV=1;
      //}else if(strcmp(argv[i], "-print_div_all")==0){
      //PRINT_DIV_ALL=1;
    }else if(strcmp(argv[i], "-print_div")==0){
      *PRINT_DIV=1;
    }else if(strcmp(argv[i], "-ali_tm")==0){
      ALI_TM=1;
    }else if(strcmp(argv[i], "-ali_co")==0){
      ALI_CO=1;
    }else if(strcmp(argv[i], "-ali_ss")==0){
      *ALI_SS=1;
    }else if(strcmp(argv[i], "-ss_mult")==0){
      *SS_MULT=1;
    }else if(strcmp(argv[i], "-shift_max")==0){
      i++; int tmp; sscanf(argv[i], "%d", &tmp);
      if(tmp>=0){*SHIFT_MAX=tmp;}
      else{printf("WARNING, %d not allowed value for SHIFT_MAX\n", tmp);}
      /*}else if(strcmp(argv[i], "-outgroup")==0){
	i++;
	if((strcmp(argv[i],"TN")!=0) && 
	(strcmp(argv[i],"CD")!=0) && 
	(strcmp(argv[i],"TM")!=0)){
	printf("WARNING, outgroup method %s not implemented, using default\n",
	argv[i]);
	}else{strcpy(OUTG, argv[i]);}
	printf("Outgroup method set to %s\n", OUTG); */
    }else if(strcmp(argv[i], "-func")==0){
      i++; sscanf(argv[i], "%s", file_fun);
    }else if(strcmp(argv[i], "-file")==0){
      i++;
    }else{
      printf("WARNING, unknown command %s", argv[i]);
    }
  }
  
  if(ali==0 && file_list[0]=='\0'){
    printf("ERROR, either MSA (-ali) or pdbs (-pdblist) must be specified\n");
    exit(8);
  }
  if(*ALI_SS){printf("Modifying MSA with sec.str.information: ");}
  if(name_in[0]=='\0'){
    if(file_list[0]){Get_file_name(name_in, file_list);}
    else if(file_ali[0]){Get_file_name(name_in, file_ali);}
    else{
      printf("ERROR, either MSA (-ali) or pdbs (-pdblist) must be specified\n");
      exit(8);
    }
  }
  if(INP_MSA){printf("Sequences input as MSA\n");}
  else{printf("Input sequences not aligned\n");}
}

float **Read_function(char *file_fun, struct Prot_input *Prot,
		      int *index, int N)
{
  if(file_fun[0]=='\0')return(NULL);
  FILE *file_in=fopen(file_fun, "r");
  if(file_in==NULL)return(NULL);

  int n=0, i, j;
  float **f_sim=malloc(N*sizeof(float *)), sim;
  for(i=0; i<N; i++){
    f_sim[i]=malloc(N*sizeof(float));
    for(j=0; j<N; j++)f_sim[i][j]=-1;
  }
  char string[1000], name1[40], name2[40];
  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(string[0]=='#')continue;
    sscanf(string, "%s%s%f", name1, name2, &sim);
    i=Find_prot(name1, Prot, index, N);
    j=Find_prot(name2, Prot, index, N);
    if((i<0)||(j<0)){
      printf("WARNING, pair %s %s not identified in function file %s\n",
	     name1, name2, file_fun);
    }else{
      f_sim[i][j]=sim; f_sim[j][i]=sim; n++;
    }
  }
  fclose(file_in);
  if(n==0){
    printf("ERROR, no pair was identified in function file %s\n",
	   file_fun);
    Empty_matrix_f(f_sim, N);
    return(NULL);
  }else{
    printf("%d protein pairs identified in %s\n", n, file_fun);
    return(f_sim);
  }
}

int Find_prot(char *name, struct Prot_input *Prot, int *index, int N){
  for(int i=0; i<N; i++){
    int i1=index[i];
    if(strncmp(name, Prot[i1].domname, 5)==0){return(i);}
  }
  return(-1);
}

float Divergence(float SI, float S0, float DMAX)
{
  float Div;
  if(SI>S0){
    Div=-log((SI-S0)/(1.-S0)); if(Div>DMAX)Div=DMAX;
  }else{
    Div=DMAX;
  }
  return(Div);
}
void Print_seq(int *ali, char *seq, int N_ali){
  for(int k=0; k<N_ali; k++){
    if(ali[k]<0){printf("-");}
    else{printf("%c",seq[ali[k]]);}
  }
  printf("\n");
}
int Count_ali(int *ali, int len){
  int n=0, k; for(k=0; k<len; k++)if(ali[k]>=0)n++;
  return(n);
}

void Summary_identical(FILE *file_id, int it_opt)
{
  char *name_c[2], *name_ct[2]; int i, j, k;
  name_c[0]=malloc(9*sizeof(char)); strcpy(name_c[0], "fct");
  name_c[1]=malloc(9*sizeof(char)); strcpy(name_c[1], "mct");
  name_ct[0]=malloc(20*sizeof(char)); strcpy(name_ct[0], "few contacts");
  name_ct[1]=malloc(20*sizeof(char)); strcpy(name_ct[1], "many contacts");
  fprintf(file_id,
	  "# ali=aligned, ide=identical aa,"
	  " sup=superimposed in space, d<d0(TM)/2"
	  " cons cont=conserved contact\n"
	  "# fct=residues with few structural contacts,"
	  " mct= many contacts\n");
  fprintf(file_id,"# %d alignments tested\n", PTYPE+1); 
  fprintf(file_id,"# npairs= %d (pdb) %d (seq)\n", npair_pdb, npair_seq);
  for(i=1; i<PTYPE; i++){
    if(i==code[2] || i==code[3] || i==code[4]) 
      fprintf(file_id,
	      "# Fraction of pairs for which %s_ali has not max %s: %.3f\n",
	      ali_name[i], ali_name[i], Diff_opt[i]/npair_pdb);
  }
  fprintf(file_id,"# Fraction of pairs where the named alignment has max PC:");
  for(i=0; i<PTYPE; i++)
    fprintf(file_id,"  %s %.3f", ali_name[i], opt_PC[i]/npair_pdb);
  fprintf(file_id,"\n");

  //
  double ave=0, se=0, re2A=0, re2B=0, ratio=0, ser=0;
  if(ALI_SS){
    fprintf(file_id,"# Secstr_corrected_alignment:\n");
    char id[3][4];
    strcpy(id[0],"SI"); strcpy(id[1],"TM"); strcpy(id[2],"CO");
    char sign[2]; sign[0]='<'; sign[1]='>';
    double *count[3];
    count[0]=sum_SS_SI; count[1]=sum_SS_TM; count[2]=sum_SS_CO;
    
    for(i=0; i<3; i++){
      fprintf(file_id,"# ");
      for(j=0; j<2; j++){
	ave=Mean_freq(&se, &re2A, count[i][j], npair_pdb);
	fprintf(file_id," P(%s%c%s(Input)= %.4f s.e.=%.4f",
		id[i],sign[j],id[i],ave, se);
      }
      fprintf(file_id,"\n");
    }
    double **count2[3]; count2[0]=sum_SS_SI_TM;
    count2[1]=sum_SS_SI_CO;count2[2]=sum_SS_TM_CO; int i2=0;
    for(i=0; i<3; i++){
      for(k=i+1; k<3; k++){
	fprintf(file_id,"# ");
	for(j=0; j<2; j++)
	  for(int j2=0; j2<2; j2++){
	    float PA=Mean_freq(&se, &re2A, count2[i2][j][j2], count[i][j]);
	    float PB=Mean_freq(&se, &re2B, count[k][j2], npair_pdb);
	    double ratio=PA/PB, ser=ratio*sqrt(re2A+re2B);
	    /*fprintf(file_id," P(%s%c,%s%c)/P(%s%c)P(%s%c)",
	      id[i],sign[j],id[k],sign[j2],
	      id[i],sign[j],id[k],sign[j2]);*/
	    fprintf(file_id," Prop(%s%c,%s%c)", id[i],sign[j],id[k],sign[j2]);
	    fprintf(file_id,"= %.2f se= %.2f",ratio, ser);
	  }
	fprintf(file_id,"\n"); i2++;
      }
    }
  } // end ALI_SS

  int npair=npair_pdb, mult=0;
  for(int it=0; it<=PTYPE; it++){
    int i=it; // i=type of alignment
    if(ALL_TYPES==0){
      if(it!=0 && it!=PTYPE){continue;}
      if(it==PTYPE){i=1;}
    }
    if(it==PTYPE){
      mult=1;
    }

    /*if(ALL_TYPES){
      if(it==PTYPE){i=it_opt; it1=it_opt; mult=1;}
      }else if(it>=CTYPE){
      break;
      }else if(it){
      i=it_opt; it1=1; mult=1;
      }*/

    fprintf(file_id,"# it= %d it_opt= %d i= %d mult= %d ALL=%d\n",
	    it, it_opt, i, mult, ALL_TYPES);

    char *ali_name_i=ali_name[it];
    if(mult){npair=npair_seq; ali_name_i=ali_name[it_opt];}

    float nind=sqrt(2*npair_seq+0.25)+0.5;
    fprintf(file_id,"######## %s alignment\n", ali_name_i);
    fprintf(file_id, "### Average pairwise scores:\n");
    fprintf(file_id, "Mean ali: %.4g\n", Sim_ave[0][i]/npair);
    fprintf(file_id, "Mean SI:  %.4f\n", Sim_ave[1][i]/npair);
    fprintf(file_id, "Mean TM:  %.4f\n", Sim_ave[2][i]/npair);
    fprintf(file_id, "Mean CO:  %.4f\n", Sim_ave[3][i]/npair);
    fprintf(file_id, "Mean PC:  %.4g\n", Sim_ave[4][i]/npair);

    fprintf(file_id, "### Differences between alignments:\n");
    for(int j=0; j<i; j++){
      if(Sim_ave[0][j]==0){continue;} double d1, d2;
      fprintf(file_id, "Diff_ali_%s_vs_%s:", ali_name[j],ali_name_i);
      d1=Ave_se(&d2, Sim_diff[1][i][j], Sim_diff_2[1][i][j], npair, nind);
      fprintf(file_id, "\tSI: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, Sim_diff[2][i][j], Sim_diff_2[2][i][j], npair, nind);
      fprintf(file_id, "\tTM: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, Sim_diff[3][i][j], Sim_diff_2[3][i][j], npair, nind);
      fprintf(file_id, "\tCO: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, Sim_diff[4][i][j], Sim_diff_2[4][i][j], npair, nind);
      fprintf(file_id, "\tPC: %.4f %.4f", d1, d2);
      d1=Ave_se(&d2, Sim_diff[0][i][j], Sim_diff_2[0][i][j], npair, nind);
      fprintf(file_id, "\tali: %.4f %.4f\n", d1, d2);
    } // end pairs of alignments


    fprintf(file_id, "### Conservation across all residues:\n");
    double sum, tot, sum_all[2], tot_all=0, re2;
    for(k=0; k<2; k++){
      sum_all[k]=Sum_bins(sum_all_ct[i][k]);
      tot_all+=sum_all[k];
    }
    if(tot_all==0){
      fprintf(file_id, "# ERROR tot_all= %.0f\n", tot_all);
      continue;
    }

    double sum_ali[2];
    for(k=0; k<2; k++)sum_ali[k]=Sum_bins(sum_ali_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_ali[0]+sum_ali[1], tot_all);
    fprintf(file_id,"Frac_ali_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(aligned)\n");  

    double sum_ide[2];
    for(k=0; k<2; k++)sum_ide[k]=Sum_bins(sum_aaid_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_ide[0]+sum_ide[1], tot_all);
    fprintf(file_id,"Frac_ide_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(identical)\n");  
    
    double sum_sup[2];
    for(k=0; k<2; k++)sum_sup[k]=Sum_bins(sum_sup_ct[i][k]);
    ave=Mean_freq(&se, &re2, sum_sup[0]+sum_sup[1], tot_all);
    fprintf(file_id,"Frac_sup_vs_all:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot_all);
    fprintf(file_id,"P(superimposed)\n");  
    
    double sum_id_cont[2], sum_ali_cont[2]; sum=0; tot=0;
    for(k=0; k<2; k++){
      sum_id_cont[k]=Sum_bins(sum_id_cont_ct[i][k]); 
      sum_ali_cont[k]=Sum_bins(sum_ali_cont_ct[i][k]);
      sum+=sum_id_cont[k]; tot+=sum_ali_cont[k];
    }
    ave=Mean_freq(&se, &re2, sum, tot);
    fprintf(file_id,"Frac_cons_cont:\t%.4f se %.4f n %.0f\t",
	    ave, se, tot);
    fprintf(file_id, "P(conserved cont) norm. by all cont\n");

    fprintf(file_id, "### Conservation many contacts vs few contacts:\n");
    // As a function of the number of contacts
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_ali[k], sum_all[k]);
      fprintf(file_id,"Frac_ali_vs_all_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_all[k]);
      fprintf(file_id,"P(aligned | %s)\n", name_ct[k]);  
    }
    
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_ide[k], sum_ali[k]);
      fprintf(file_id,"Frac_ide_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(identical | aligned, %s)\n", name_ct[k]);  
    }

    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_sup[k], sum_ali[k]);
      fprintf(file_id,"Frac_sup_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(superimposed | aligned, %s)\n", name_ct[k]);  
    }
 
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_id_cont[k], sum_ali_cont[k]);
      fprintf(file_id,
	      "Frac_cons_cont_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali_cont[k]);
      fprintf(file_id,"P(conserved cont | aligned cont, %s)\n", name_ct[k]);  
    }       

    // Neighbors in space
    fprintf(file_id, "### Spatial neighbors even if not aligned:\n");
    double sum_neigh_ali[2], sum_neigh_noali[2], tot_neigh[2]; sum=0;
    for(k=0; k<2; k++){
      sum_neigh_ali[k]=Sum_bins(sum_neigh_ali_ct[i][k]);
      sum_neigh_noali[k]=Sum_bins(sum_neigh_noali_ct[i][k]);
      tot_neigh[k]=sum_neigh_ali[k]+sum_neigh_noali[k];
      ave=Mean_freq(&se, &re2, tot_neigh[k], sum_ali[k]);
      fprintf(file_id,"Frac_neigh_vs_ali_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_ali[k]);
      fprintf(file_id,"P(neigh| %s) aligned or not\n", name_ct[k]);  
    }
    
    for(k=0; k<2; k++){
      ave=Mean_freq(&se, &re2, sum_neigh_noali[k], tot_neigh[k]);
      fprintf(file_id,"Frac_neigh_noali_vs_neigh_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, tot_neigh[k]);
      fprintf(file_id,"P(neigh, not ali| neigh, %s)\n", name_ct[k]);
      /*sum=0; for(j=0; j<=IBIN; j++)sum+=sum_neigh_noali_aaid_ct[i][k][j];
      ave=Mean_freq(&se, &re2, sum, sum_neigh_noali[k]);
      fprintf(file_id,
	    "Frac_neigh_noali_id_ali_vs_neigh_noali:\t%.4f se %.4f n %.0f\t",
	    ave, se, sum);
	    fprintf(file_id,"P(id ali | neighb not ali)\n");  */
    }

    for(k=0; k<2; k++){
      sum=Sum_bins(sum_sup_noneigh_ct[i][k]);
      ave=Mean_freq(&se, &re2, sum, sum_sup[k]);
      fprintf(file_id,"Frac_sup_noneigh_vs_sup_%s:\t%.4f se %.4f n %.0f\t",
	      name_c[k], ave, se, sum_sup[k]);
      double sum_shift=Sum_bins(sum_shift_sup_noneigh_ct[i][k]);
      if(sum)sum_shift/=sum;
      fprintf(file_id,"shift=%.1f\t",sum_shift);
      fprintf(file_id,"P(no neigh | ali sup, %s)\n", name_ct[k]);
    }

    // Propensity superimposed-identical
    fprintf(file_id, "### Propensity superimposed - identical:\n");
    double sum_id_sup[2], sum_noid_nosup[2];
    for(k=0; k<2; k++){
      sum_id_sup[k]=Sum_bins(sum_aaid_sup_ct[i][k]);
      sum_noid_nosup[k]=Sum_bins(sum_noid_nosup_ct[i][k]);
    }

    for(k=0; k<2; k++){
      double P_sup_id=Mean_freq(&se, &re2A, sum_id_sup[k], sum_ide[k]);
      double sup_noid=sum_sup[k]-sum_id_sup[k];   
      double noid=sum_ali[k]-sum_ide[k];
      double P_sup_noid=Mean_freq(&se, &re2B, sup_noid, noid);
      ratio=P_sup_id/P_sup_noid; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"P(sup|ide)/P(sup|no id)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, sum_ide[k], noid);

      double P_nosup_noid=Mean_freq(&se, &re2A, sum_noid_nosup[k], noid);
      double nosup=sum_ali[k]-sum_sup[k];
      double nosup_id=nosup-sum_noid_nosup[k];   
      double P_nosup_id=Mean_freq(&se, &re2B, nosup_id, sum_ide[k]);
      ratio=P_nosup_noid/P_nosup_id; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"P(no sup|no id)/P(no sup|id)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, noid, sum_ide[k]);

      double P_id_sup=Mean_freq(&se, &re2A, sum_id_sup[k], sum_sup[k]);
      double id_nosup=sum_ide[k]-sum_id_sup[k];
      double P_id_nosup=Mean_freq(&se, &re2B, id_nosup, nosup);
      ratio=P_id_sup/P_id_nosup; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"P(ide|sup)/P(ide|no sup)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, sum_sup[k], nosup);

      double P_noid_nosup=Mean_freq(&se, &re2A, sum_noid_nosup[k], nosup);
      double noid_sup=noid-sum_noid_nosup[k];
      double P_noid_sup=Mean_freq(&se, &re2B, noid_sup, sum_sup[k]);
      ratio=P_noid_nosup/P_noid_sup; ser=ratio*sqrt(re2A+re2B);
      fprintf(file_id,"P(no id|no sup)/P(no id|sup)_%s:", name_c[k]);
      fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, nosup, sum_sup[k]);

    } 

    // Conserved contacts
    double id_cont_sup[3], ali_cont_sup[3];
    double id_cont_aaid[3], ali_cont_aaid[3];
    double id_cont=0, ali_cont=0;
    for(k=0; k<3; k++){
      id_cont_sup[k]=0; ali_cont_sup[k]=0; 
      id_cont_aaid[k]=0; ali_cont_aaid[k]=0; 
      for(j=0; j<=IBIN; j++){
	id_cont_sup[k]+=sum_id_cont_sup[i][k][j];
	ali_cont_sup[k]+=sum_ali_cont_sup[i][k][j];
	id_cont_aaid[k]+=sum_id_cont_aaid[i][k][j];
	ali_cont_aaid[k]+=sum_ali_cont_aaid[i][k][j];
      }
      id_cont+=id_cont_sup[k];
      ali_cont+=ali_cont_sup[k];
    }
    // Conditional Prob conserved contact | superimposed
    fprintf(file_id, "### Propensity ide or sup | ide or sup in contact\n");

    double sum_a=sum_ali[0]+sum_ali[1]; //0=few cont 1=many cont
    double sum_i=sum_ide[0]+sum_ide[1];
    double sum_s=sum_sup[0]+sum_sup[1];
    double re2i,Pi=Mean_freq(&se, &re2i, sum_i, sum_a);//re2noi=re2i*(1-Pi)/Pi;
    double re2s,Ps=Mean_freq(&se, &re2s, sum_s, sum_a);//re2nos=re2s*(1-Ps)/Ps;
    
    sum_i=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    double sum_a2=2*(ali_cont_aaid[0]+ali_cont_aaid[1]+ali_cont_aaid[2]);
    double re2ic=0, Pic=Mean_freq(&se, &re2ic, sum_i, sum_a2),
      re2noic=re2ic*(1-Pic)/Pic;
    fprintf(file_id, "P(id|cont)/P(id)= %.4f se %.2g\n", Pic/Pi, se/Pi);
    sum_s=ali_cont_sup[1]+2*ali_cont_sup[2];
    sum_a2=2*(ali_cont_sup[0]+ali_cont_sup[1]+ali_cont_sup[2]);
    double re2sc=0, Psc=Mean_freq(&se, &re2sc, sum_s, sum_a2),
      re2nosc=re2sc*(1-Psc)/Psc;
    fprintf(file_id, "P(sup|cont)/P(sup)= %.4f se %.2g\n", Psc/Ps, se/Ps);

    // Correlated substitutions
    sum=2*ali_cont_aaid[0]+ali_cont_aaid[1];
    double PA=Mean_freq(&se, &re2A, 2*ali_cont_aaid[0], sum);
    ratio=PA/(1-Pic); ser=ratio*sqrt(re2A+re2noic);
    fprintf(file_id,"P(no id|cont no id)/P(no id):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n",ratio,ser, sum, sum_a);

    //sum=2*ali_cont_aaid[0]+ali_cont_aaid[1];
    PA=Mean_freq(&se, &re2A, ali_cont_aaid[1], sum);
    ratio=PA/Pic; ser=ratio*sqrt(re2A+re2ic);
    fprintf(file_id,"P(id|cont no id)/P(id):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_aaid[2], sum);
    ratio=PA/Pic; ser=ratio*sqrt(re2A+re2ic);
    fprintf(file_id,"P(id|cont id)/P(id):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=ali_cont_aaid[1]+2*ali_cont_aaid[2];
    PA=Mean_freq(&se, &re2A, ali_cont_aaid[1], sum);
    ratio=PA/(1-Pic); ser=ratio*sqrt(re2A+re2noic);
    fprintf(file_id,"P(no id|cont id)/P(no id):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=2*ali_cont_sup[0]+ali_cont_sup[1];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_sup[0], sum);
    ratio=PA/(1-Psc); ser=ratio*sqrt(re2A+re2nosc);
    fprintf(file_id,"P(no sup|cont no sup)/P(no sup):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=2*ali_cont_sup[0]+ali_cont_sup[1];
    PA=Mean_freq(&se, &re2A, ali_cont_sup[1], sum);
    ratio=PA/Psc; ser=ratio*sqrt(re2A+re2sc);
    fprintf(file_id,"P(sup|cont no sup)/P(sup):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    sum=ali_cont_sup[1]+2*ali_cont_sup[2];
    PA=Mean_freq(&se, &re2A, 2*ali_cont_sup[2], sum);
    ratio=PA/Psc; ser=ratio*sqrt(re2A+re2sc);
    fprintf(file_id,"P(sup|cont sup)/P(sup):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    //sum=ali_cont_sup[1]+2*ali_cont_sup[2];
    PA=Mean_freq(&se, &re2A, ali_cont_sup[1], sum);
    ratio=PA/(1-Psc); ser=ratio*sqrt(re2A+re2nosc);
    fprintf(file_id,"P(no sup|cont sup)/P(no sup):");
    fprintf(file_id,"\t%.4f se %.2g n %.0f %.0f\n", ratio,ser, sum, sum_a);

    fprintf(file_id, "### Conditional probability conserved contact | sup:\n");
    double P_cons_cont_sup[3], re2ccs[3];
    for(k=0; k<3; k++){
      P_cons_cont_sup[k]=
	Mean_freq(&se, re2ccs+k, id_cont_sup[k], ali_cont_sup[k]);
      fprintf(file_id,"Frac_cons_cont_vs_%d_sup:\t%.4f se %.4f n %.0f\t",
	      k, P_cons_cont_sup[k], se, ali_cont_sup[k]);
      fprintf(file_id, "P(conserved cont | %d sup res)\n",k);
    }

    // Conditional Prob conserved contact | identical 
    fprintf(file_id,
	    "### Conditional probability conserved contact | id aa:\n");
    double P_cons_cont_aaid[3], re2cca[3];
    for(k=0; k<3; k++){
      P_cons_cont_aaid[k]=
	Mean_freq(&se, re2cca+k, id_cont_aaid[k], ali_cont_aaid[k]);
      fprintf(file_id,"Frac_cons_cont_vs_%d_id:\t%.4f se %.4f n %.0f\t",
	      k, P_cons_cont_aaid[k], se, ali_cont_aaid[k]);
      fprintf(file_id, "P(conserved cont | %d identical aa)\n",k);
    }


    // Propensity conserved contacts - identical
    fprintf(file_id, "### Propensity conserved contact - identical:\n");
    ratio=P_cons_cont_aaid[2]/(P_cons_cont_aaid[0]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"P(cons cont|2 ide)/P(cons cont|0 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[2], ali_cont_aaid[0]+ali_cont_aaid[1]);
    //
    //
    ratio=(1-P_cons_cont_aaid[1])/(1-P_cons_cont_aaid[2]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"P(chg cont|1 ide)/P(chg cont|2 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[1], ali_cont_aaid[2]);
    //
    ratio=(1-P_cons_cont_aaid[0])/(1-P_cons_cont_aaid[2]);
    ser=ratio*sqrt(re2cca[2]+re2cca[0]);
    fprintf(file_id,"P(chg cont|0 ide)/P(chg cont|2 ide):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_aaid[0], ali_cont_aaid[2]);
    //
    double P_2aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[2], id_cont);
    double P_2aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[2]-id_cont_aaid[2], ali_cont-id_cont);
    ratio=P_2aaid_id_cont/P_2aaid_ch_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(2 ide|cons cont)/P(2 ide|chg cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, id_cont, ali_cont-id_cont);
    //
    double P_1aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[1], id_cont);
    double P_1aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[1]-id_cont_aaid[1], ali_cont-id_cont);
    ratio=P_1aaid_ch_cont/P_1aaid_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(1 ide|chg cont)/P(1 ide|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
    ratio, ser, ali_cont-id_cont, id_cont);
    //
    double P_0aaid_id_cont=Mean_freq(&se, &re2A, id_cont_aaid[0], id_cont);
    double P_0aaid_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_aaid[0]-id_cont_aaid[0], ali_cont-id_cont);
    ratio=P_0aaid_ch_cont/P_0aaid_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(0 ide|chg cont)/P(0 ide|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);

    // Propensity conserved contacts - superimposed
    fprintf(file_id, "### Propensity conserved contact - superimposed:\n");
    ratio=P_cons_cont_sup[2]/P_cons_cont_sup[0];
    ser=ratio*sqrt(re2ccs[2]+re2ccs[0]);
    fprintf(file_id,"P(cons cont|2 sup)/P(cons cont|0 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[2], ali_cont_sup[0]);
    //
    ratio=(1-P_cons_cont_sup[1])/(1-P_cons_cont_sup[2]);
    ser=ratio*sqrt(re2ccs[2]+re2ccs[1]);
    fprintf(file_id,"P(chg cont|1 sup)/P(chg cont|2 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[1], ali_cont_sup[2]);
    //
    ratio=(1-P_cons_cont_sup[0])/(1-P_cons_cont_sup[2]);
    ser=ratio*sqrt(re2ccs[2]+re2ccs[0]);
    fprintf(file_id,"P(chg cont|0 sup)/P(chg cont|2 sup):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont_sup[0], ali_cont_sup[2]);
    //
    double P_2sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[2], id_cont);
    double P_2sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[2]-id_cont_sup[2], ali_cont-id_cont);
    ratio=P_2sup_id_cont/P_2sup_ch_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(2 sup|cons cont)/P(2 sup|chg cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, id_cont, ali_cont-id_cont);
    //
    double P_1sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[1], id_cont);
    double P_1sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[1]-id_cont_sup[1], ali_cont-id_cont);
    ratio=P_1sup_ch_cont/P_1sup_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(1 sup|chg cont)/P(1 sup|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);
    //
    double P_0sup_id_cont=Mean_freq(&se, &re2A, id_cont_sup[0], id_cont);
    double P_0sup_ch_cont=
      Mean_freq(&se, &re2B,ali_cont_sup[0]-id_cont_sup[0], ali_cont-id_cont);
    ratio=P_0sup_ch_cont/P_0sup_id_cont; ser=ratio*sqrt(re2A+re2B);
    fprintf(file_id,"P(0 sup|chg cont)/P(0 sup|cons cont):");
    fprintf(file_id,"\t%.4f se %.4f n %.0f %.0f\n",
	      ratio, ser, ali_cont-id_cont, id_cont);

    fprintf(file_id, "### As a function of sequence identity:\n");
    fprintf(file_id, "#1=f_id s.e. 3=f_ali s.e. ");
    fprintf(file_id, " 5=f_sup s.e. 7=f_cons_cont s.e.");
    fprintf(file_id, " 9=f_id_mct/f_id_fct s.e.");
    fprintf(file_id, " 11=f_sup_mct/f_sup_fct s.e.");
    fprintf(file_id, " 13=f_cons_cont_mct/f_cons_cont_fct s.e.");
    fprintf(file_id, " 15=f_neigh_noali_mct/fct s.e.");
    fprintf(file_id, " 17=f_id_sup/f_id*f_sup s.e.");
    fprintf(file_id, " 19=f_noid_nosup/f_noid*f_nosup s.e.");
    fprintf(file_id, " 21=f_id_cons_cont/f_id*f_cons_cont s.e.");
    fprintf(file_id, " 23=f_noid_nocons_cont/f_noid*f_nocons_cont s.e.");
    fprintf(file_id, " 25=f_sup_cons_cont/f_sup*f_cons_cont s.e.");
    fprintf(file_id, " 27=f_nosup_nocons_cont/f_nosup*f_nocons_cont s.e.");
    /*for(int k=0; k<3; k++){
      fprintf(file_id, " %d=f_cont_%dsup s.e.", 23+2*k, k);
      }*/
    fprintf(file_id, " 29=n\n");

    float P_A, P_B;
    for(j=0; j<=IBIN; j++){
      double ali=sum_ali_ct[i][0][j]+sum_ali_ct[i][1][j];
      double all=sum_all_ct[i][0][j]+sum_all_ct[i][1][j];
      if(ali<10)continue;
      double aaid=sum_aaid_ct[i][0][j]+sum_aaid_ct[i][1][j];
      double re2id=0, P_aaid=Mean_freq(&se, &re2id, aaid, ali);
      fprintf(file_id, "%.3f\t%.1g", P_aaid, se);
      ave=Mean_freq(&se, &re2, ali, all);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      double sup=sum_sup_ct[i][0][j]+sum_sup_ct[i][1][j];
      double re2sup, P_sup=Mean_freq(&se, &re2sup, sup, ali);
      fprintf(file_id, "\t%.3f\t%.1g", P_sup, se);
      double cons_cont=sum_id_cont_ct[i][0][j]+sum_id_cont_ct[i][1][j];
      double ali_cont=sum_ali_cont_ct[i][0][j]+sum_ali_cont_ct[i][1][j];
      double re2cont=0, P_cont=Mean_freq(&se, &re2cont, cons_cont, ali_cont);
      fprintf(file_id, "\t%.3f\t%.1g", P_cont, se);

      P_A=Mean_freq(&se, &re2A, sum_aaid_ct[i][1][j], sum_ali_ct[i][1][j]);
      P_B=Mean_freq(&se, &re2B, sum_aaid_ct[i][0][j], sum_ali_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);     
      P_A=Mean_freq(&se, &re2A, sum_sup_ct[i][1][j], sum_ali_ct[i][1][j]);
      P_B=Mean_freq(&se, &re2B, sum_sup_ct[i][0][j], sum_ali_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      P_A=Mean_freq(&se,&re2A,sum_id_cont_ct[i][1][j],sum_ali_cont_ct[i][1][j]);
      P_B=Mean_freq(&se,&re2B,sum_id_cont_ct[i][0][j],sum_ali_cont_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      // Neighbors not aligned, mct/fct
      P_A=Mean_freq(&se,&re2A,sum_neigh_noali_ct[i][1][j],sum_all_ct[i][1][j]);
      P_B=Mean_freq(&se,&re2B,sum_neigh_noali_ct[i][0][j],sum_all_ct[i][0][j]);
      ratio=P_A/P_B; se=ratio*sqrt(re2A+re2B);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);

      // Propensities id-sup
      double id_sup=sum_aaid_sup_ct[i][0][j]+sum_aaid_sup_ct[i][1][j];
      double P_id_sup=Mean_freq(&se, &re2A, id_sup, sup);
      ratio=P_id_sup/P_aaid; se=ratio*sqrt(re2A+re2id);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      double noid_nosup=sum_noid_nosup_ct[i][0][j]+sum_noid_nosup_ct[i][1][j];
      double P_noid_nosup=Mean_freq(&se, &re2A, noid_nosup, ali-sup);
      ratio=P_noid_nosup/(1-P_aaid); se=ratio*sqrt(re2A+re2id);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      // Propensities id-cons cont
      ratio=P_cons_cont_aaid[2]/P_cont; se=ratio*sqrt(re2cont+re2cca[2]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      ratio=(1-P_cons_cont_aaid[0])/(1-P_cont);
      se=ratio*sqrt(re2cont+re2cca[0]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      // Propensities sup-cons cont
      ratio=P_cons_cont_sup[2]/P_cont; se=ratio*sqrt(re2cont+re2ccs[2]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);
      ratio=(1-P_cons_cont_sup[0])/(1-P_cont);
      se=ratio*sqrt(re2cont+re2ccs[0]);
      fprintf(file_id, "\t%.3f\t%.1g", ratio, se);

      /*float sum_aaid_sup=sum_aaid_sup_ct[i][0][j]+sum_aaid_sup_ct[i][1][j];
      ave=Mean_freq(&se, &re2, sum_aaid_sup, sum_aaid);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      float sup_noaa=sum_sup_ct[i][0][j]+sum_sup_ct[i][1][j]-sum_aaid_sup;
      float noaa=sum_ali-sum_aaid;
      ave=Mean_freq(&se, &re2, sup_noaa, noaa);
      fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      for(k=0; k<3; k++){
	ave=Mean_freq(&se,&re2,sum_id_cont_sup[i][k][j],
		      sum_ali_cont_sup[i][k][j]);
	fprintf(file_id, "\t%.3f\t%.1g", ave, se);
      }*/
      fprintf(file_id, "\t%.0f\n", ali);
    }
      
  }	    

  for(int k=0; k<2; k++){free(name_c[k]); free(name_ct[k]);}
}

void Count_noali(int *noali, int *ali,  int *ncont, int n, float c_ave){
  noali[0]=0; noali[1]=0;
  for(int i=0; i<n; i++){
    if(ali[i]<0){if(ncont[i]<=c_ave){noali[0]++;}else{noali[1]++;}}
  }
}
double Sum_bins(double *bin){
  double sum=0; for(int j=0; j<=IBIN; j++)sum+=bin[j];
  return(sum);
}

float Ave_se(double *se, double sum1, double sum2, int n, float nind){
  double ave=sum1/n;
  *se=sqrt((sum2/n-ave*ave)/nind);
  return(ave);
}


void Sec_str_AA_propensity(struct protein *pi, struct protein *pj,
			   int *ali, int nali, float DP)
{
  int ns=4, na=22, gap=3, gap2=21, nbin=5, nmax=nbin-1;
  float lbin=0.25;
  if(Prop_secstr==NULL){
    Prop_secstr=malloc(nbin*sizeof(double **));
    Prop_AA=malloc(nbin*sizeof(double **));
    DP_bin=malloc(nbin*sizeof(double));
    DP_norm=malloc(nbin*sizeof(double));
    for(int i=0; i<nbin; i++){
      Prop_secstr[i]=Allocate_mat2_d(ns, ns);
      Prop_AA[i]=Allocate_mat2_d(na, na);
      DP_bin[i]=0; DP_norm[i]=0;
    }
    for(int i=0; i<4; i++)SS_name[i]=malloc(10*sizeof(char));
    strcpy(SS_name[0], "Coil");
    strcpy(SS_name[1], "Helix");
    strcpy(SS_name[2], "Strand");
    strcpy(SS_name[3], "Gap");
    for(int i=0; i<20; i++)AA_name[i]=AANAME1[i];
    AA_name[20]='-'; AA_name[21]='X'; 
  }

  int bin=DP/lbin; if(bin>nmax)bin=nmax;
  DP_bin[bin]+=DP; DP_norm[bin]++;
  double **Prop_bin=Prop_secstr[bin];
  double **Prop_AA_bin=Prop_AA[bin];

  int i, j_old=0;
  for(i=0; i<pi->len; i++){
    int iss1=pi->ss3[i];
    int aa1=Code_AA_2(pi->aseq[i], AANAME1, 20);
    int j=ali[i];

    if(j<0){
      Prop_bin[iss1][gap]++;
      Prop_AA_bin[aa1][gap2]++;
    }else{
      int j1=j-1;
      if(j_old<j1){
	for(int jj=j_old; jj<j1; jj++){
	  int iss2=pj->ss3[jj];
	  int aa2=Code_AA_2(pj->aseq[jj], AANAME1, 20);
	  Prop_bin[iss2][gap]++;
	  Prop_AA_bin[aa2][gap2]++;
	}
      }
      j_old=j;
      int iss2=pj->ss3[j];
      int aa2=Code_AA_2(pj->aseq[j], AANAME1, 20);
      Prop_bin[iss1][iss2]++;
      Prop_AA_bin[aa1][aa2]++;
    }
  }
}

void Print_propensities(char *nameout)
{
  int ns=4, nbin=5; // na=22, gap=3, gap2=21, nmax=nbin-1;
  //float lbin=0.25;
  FILE *file_out=fopen(nameout, "w");
  printf("Writing Sec.str. and AA propensities of aligned sites in %s\n",
	 nameout);
  char name_ss[200], tmp[20]; int a;
  sprintf(name_ss, "%s", SS_name[0]);
  for(a=1; a<ns; a++){sprintf(tmp, "\t%s", SS_name[a]); strcat(name_ss, tmp);}
  char name_AA[200];
  sprintf(name_AA, "%c", AA_name[0]);
  for(a=1; a<20; a++){sprintf(tmp, "\t%c", AA_name[a]); strcat(name_AA, tmp);}
  strcat(name_AA, "\t-\n");

  for(int i=0; i<nbin; i++){
    fprintf(file_out, "# Bin %d: <PC_div>= %.3g norm= %.3g\n",
	    i+1, DP_bin[i]/DP_norm[i], DP_norm[i]);
    if(DP_norm[i]<4)continue;
    int n=ns, a, b; double **Prop=Prop_secstr[i], F[n];
    // Symmetrize
    for(a=0; a<n; a++){ 
      for(b=0; b<a; b++){
	Prop[a][b]+=Prop[b][a];
	Prop[a][b]/=2;
	Prop[b][a]=Prop[a][b];
      }
    }
    int zero=0; double norm=0;
    for(a=0; a<n; a++){
      F[a]=0; for(b=0; b<n; b++)F[a]+=Prop[a][b];
      norm+=F[a]; if(F[a]==0)zero=1;
    }

    float Low=-log(10);
    fprintf(file_out, "# Secondary structure: Frequency\n");
    fprintf(file_out, "# %s\n", name_ss);
    fprintf(file_out, "# %.3f", F[0]/norm);
    for(a=1; a<n; a++)fprintf(file_out, "\t%.3f", F[a]/norm);
    fprintf(file_out, "\n");
    if(zero)continue;
    fprintf(file_out, "# Secondary structure: Propensities\n");
    fprintf(file_out, "#SS\t%s\n", name_ss);
    for(a=0; a<n; a++){
      fprintf(file_out, "%s", SS_name[a]);
      for(b=0; b<n; b++){
	double p=norm*Prop[a][b]/(F[a]*F[b]);
	if(p>0){p=log(p);}else{p=Low;}
	fprintf(file_out, "\t%.4g", p);
      }
      fprintf(file_out, "\n");
    }
  }
  fclose(file_out);
}

char **Assign_Seq(char ***name_seq, int *len_seq,
		  struct Prot_input *Prot_in, int N_seq,
		  int *rep_str, int *i_seq, int N_ali)
{
  *name_seq=malloc(N_seq*sizeof(char *));
  char **Seq=malloc(N_seq*sizeof(char *));
  for(int i=0; i<N_seq; i++){
    struct Prot_input *proti=Prot_in+i_seq[rep_str[i]];
    (*name_seq)[i]=malloc(200*sizeof(char));
    //strcpy((*name_seq)[i], proti->domname);
    sprintf((*name_seq)[i], "%s %c %s",
	    proti->domname, proti->chain, proti->domain);
    int Li=proti->len, k=0;
    //printf("Assigning sequence %d L= %d\n", i, Li);
    len_seq[i]=Li;
    Seq[i]=malloc(Li* sizeof(char));
    for(int j=0; j<N_ali; j++){
      if(proti->seq[j]!='-' && k<Li){Seq[i][k]=proti->seq[j]; k++;}
    }
    if(k!=Li){
      printf("ERROR in Assign_seq, L= %d %d\n", Li, k); exit(8);
    }
  }
  return(Seq);
}

void Write_ali_pair(int ***Ali_pair, int i, int j, int *L_seq, int *ali_PC)
{
  int *ali_ij=Ali_pair[i][j], *ali_ji=NULL, u;
  int Li=L_seq[i], Lj=L_seq[j];
  if(ALL_PAIRS){
    ali_ji=Ali_pair[j][i]; for(u=0; u<Lj; u++){ali_ji[u]=-1;}
  }
  int error=0, *ali=ali_PC;
  for(u=0; u<Li; u++){
    ali_ij[u]=*ali;
    if(*ali>=Lj){
      ali_ij[u]=-1;
      if(error)continue;
      printf("ERROR in write_ali_pair %d %d, res %d >= %d\n",
	     i, j, *ali, Lj);
      for(int v=0; v<Li; v++)printf("%d ",ali_PC[v]);
      printf("\n"); error=1;
    }else if(ali_ji && *ali>=0){
      ali_ji[*ali]=u;
    }
    ali++;
  }
}

int ***Select_alis(int ***Ali_pair, int N_seq, int *len_seq, int *rep_str)
{
  int ***Ali_pair_seq=malloc(N_seq*sizeof(int **));
  int n=N_seq;
  for(int i=0; i<N_seq; i++){
    if(ALL_PAIRS==0){n=i; if(n==0)continue;}
    Ali_pair_seq[i]=malloc(n*sizeof(int *));
    int Li=len_seq[i], i_str=rep_str[i];
    for(int j=0; j<n; j++){
      Ali_pair_seq[i][j]=malloc(Li*sizeof(int));
      if(j==i)continue;
      int *ali_sel=Ali_pair_seq[i][j], j_str=rep_str[j];
      int i1, i2;
      if(i_str<j_str){i1=i_str; i2=j_str;}
      else{i1=j_str; i2=i_str;}
      int *ali=Ali_pair[i2][i1], k;
      if(i2==i_str){
	for(k=0; k<Li; k++)ali_sel[k]=ali[k];
      }else{
	for(k=0; k<Li; k++)ali_sel[k]=-1;
	for(k=0; k<len_seq[j]; k++)
	  if(ali[k]>=0 && ali[k]<Li)ali_sel[ali[k]]=k;
      }
    }
  }
  return(Ali_pair_seq);
}

int *Representative_structure(int *N_conf, int **conformation, int n,
			      struct protein *prots, float **PC_all)
{
  int *rep_str=malloc(n*sizeof(int));
  for(int i=0; i<n; i++){
    if(N_conf[i]==1){
      rep_str[i]=conformation[i][0];
    }else if(N_conf[i]==2){
      int c1=conformation[i][0], c2=conformation[i][1];
      if(prots[c1].len > prots[c2].len){rep_str[i]=c1;}else{rep_str[i]=c2;}
    }else{
      int nc=N_conf[i];
      float sim_conf_ave[nc];
      for(int j=0; j<nc; j++){
	sim_conf_ave[j]=0;
	int cj=conformation[i][j];
	for(int k=0; k<j; k++){
	  int ck=conformation[i][k], c1, c2;
	  if(cj>ck){c2=cj; c1=ck;}else{c2=ck; c1=cj;}
	  sim_conf_ave[j]+=PC_all[c2][c1];
	  sim_conf_ave[k]+=PC_all[c2][c1];
	}
      }
      int jmax=0;
      for(int j=1; j<nc; j++)
	if(sim_conf_ave[j]>sim_conf_ave[jmax])jmax=j;
      rep_str[i]=conformation[i][jmax];
    }
  }

  return(rep_str);
}

void Set_type(char **what, int *code, int *ATYPE, char *name, int it){
  what[*ATYPE]=malloc(10*sizeof(char));
  strcpy(what[*ATYPE], name);
  code[it]=*ATYPE;
  (*ATYPE)++;
}

FILE *Open_file_div(char *name_div, char *name_in, char *head)
{
  Change_ext(name_div, name_in, ".prot.div");
  FILE *file_div=fopen(name_div, "w");
  printf("Printing grouped structure divergence in %s\n", name_div);
  fprintf(file_div, "%s", head);
  fprintf(file_div, "###Prot1 Prot2 "); int k=3;
  if(INP_MSA){
    fprintf(file_div,
	    " %d=TN_Div_INP %d=Cont_Div_INP %d=TM_Div_INP %d=PC_Div_INP",
	    k, k+1, k+2, k+3); k+=4;
  }
  fprintf(file_div," %d=TN_Div_PC %d=Cont_Div_PC %d=TM_Div_PC %d=PC_Div_PC\n",
	  k, k+1, k+2, k+3);
  fprintf(file_div, "### 0 0");
  if(INP_MSA)fprintf(file_div, "  1 1 1 1");
  fprintf(file_div, "  1 1 1 1\n");
  return(file_div);
}

FILE *Open_file_sim(char *name_sim, char *name_in, char *head, float **fun_sim)
{
  Change_ext(name_sim, name_in, ".prot.sim");
  FILE *file_sim=fopen(name_sim, "w");
  printf("Printing grouped structure similarity in %s\n", name_sim);

  fprintf(file_sim, "%s", head);
  fprintf(file_sim, "### Prot1 Prot2 "); int k=3;
  if(INP_MSA){
    fprintf(file_sim,
	    " %d=nali_INP %d=SI_INP %d=CO_INP %d=TM_INP %d=PC_INP",
	    k, k+1, k+2, k+3, k+4); k+=5;
  }
  fprintf(file_sim,
	  " %d=nali_PC %d=SI_PC %d=CO_PC %d=TM_PC %d=PC_PC",
	  k, k+1, k+2, k+3, k+4);
  fprintf(file_sim, " %d=RMSD_PC", k+5);
  if(fun_sim)fprintf(file_sim, " %d=Function_similarity", k+6);
  fprintf(file_sim,"\n");
  fprintf(file_sim, "### 0 0");
  if(INP_MSA)fprintf(file_sim, "  1 1 1 1 0");
  fprintf(file_sim, "  1 1 1 1 1");
  fprintf(file_sim, " 1");
  if(fun_sim)fprintf(file_sim, " 0");
  fprintf(file_sim,"\n");
  return(file_sim);
}

float Normalization(float *norm_c, struct protein *proti, struct protein *protj)
{
  int norm_ali;
  if(NORM==0){ // Smaller length
    if(protj->len<proti->len){norm_ali=protj->len;}
    else{norm_ali=proti->len;}
    if(protj->N_cont<proti->N_cont){*norm_c=protj->N_cont;}
    else{*norm_c=proti->N_cont;}
  }else if(NORM==1){ // longer length
    if(protj->len>proti->len){norm_ali=protj->len;}
    else{norm_ali=proti->len;}
    if(protj->N_cont>proti->N_cont){*norm_c=protj->N_cont;}
    else{*norm_c=proti->N_cont;}
  }else{ // Geometric average
    norm_ali=sqrt(proti->len*protj->len);
    *norm_c=sqrt(proti->N_cont*protj->N_cont);
  }
  return(norm_ali);
}
 
 void Score_alignment(float *nali, float *SI, float *TM, float *CO,
		     float *PC, int i, int j, float **Seq_diff,
		     float **na_all, float **SI_all, float **TM_all,
		     float **CO_all, float **PC_all,
		     float **Rot, float *Shift,
		     float c_ave, float *d02, float **d2,
		     int *ali_all, int *ali_ij, int Comp_TM, float TMs,
		     struct protein *proti, struct protein *protj,
		     float norm_c, float norm_ali, int it, int last)
{
  int comp_id=0;
  if(PRINT_ID && ((ALL_TYPES && it<PTYPE) || it==0 || last)){comp_id=1;}
  
  if(ali_all && ali_all!=ali_ij){
    for(int s=0; s<proti->len; s++)ali_all[s]=ali_ij[s];
  }

  nali[it]=Count_ali(ali_ij, proti->len)/norm_ali;

  int diff=0;
  SI[it]=Seqid(&diff, ali_ij, id_aa[it],
	       proti->aseq, proti->len,
	       protj->aseq, protj->len);
  if(diff<Seq_diff[i][j]){Seq_diff[i][j]=diff;}

  // TM score and statistics of optimal rotation
  if(Comp_TM){
    TM[it]=TM_score(d2, d02, Rot, Shift, ali_ij, norm_ali, proti, protj, 0);
    //if(NORMA){TM[it]/=nali[it];}
    if(it==0 && INP_MSA && TMs>0)TM[it]=TM_all[i][j];
  }else if(d2 && comp_id){
      All_distances(d2,proti->xca_rot,proti->len,protj->xca_rot,protj->len);
  }
  CO[it]=Contact_overlap(ali_ij,
			 proti->Cont_map, proti->len,
			 protj->Cont_map, protj->len);

  PC[it]=(PC_load[0]*nali[it]+PC_load[1]*SI[it]+
	  PC_load[2]*TM[it]+PC_load[3]*CO[it])/PC_norm;

  if(na_all){
    na_all[i][j]=nali[it];
    SI_all[i][j]=SI[it];
    TM_all[i][j]=TM[it];
    CO_all[i][j]=CO[it];
    PC_all[i][j]=PC[it];
  }

  if(comp_id){
    int mult=0, it1=it;
    if(it>=PTYPE){
      mult=1; if(ALL_TYPES){it1=PTYPE;}else{it1=1;}
    }

    /*// Compare alignments
    int c1, c2; // ii, jj;
    float ***sim2[5]; sim2[0]=na_all; sim2[1]=SI_all;
    sim2[2]=TM_all; sim2[3]=CO_all; sim2[4]=PC_all;

    for(int jt=0; jt<it1; jt++){
      int ii, jj;
      if(mult && jt<PTYPE){ii=c2; jj=c1;}else{ii=i; jj=j;}
      for(int k=0; k<5; k++){
	float d=sim[k][it]-sim2[k][jt][ii][jj];
	Sim_diff[k][it1][jt]+=d;
	Sim_diff_2[k][it1][jt]+=d*d;
      }
      }*/

    if(d2){
      Examine_neighbors(d2, ali_ij, *d02, shift[it1], id_sup[it1],
			neigh_ali[it1],neigh_noali[it1],neigh_noali_aaid[it1],
			proti->aseq, proti->len, protj->aseq, protj->len);
    }

    Test_contacts(ali_ij, id_sup[it1], c_ave,
		  ali_cont_ct[it1],  id_cont_ct[it1],
		  ali_cont_sup[it1], id_cont_sup[it1],
		  ali_cont_aaid[it1],id_cont_aaid[it1],
		  proti->Cont_map, proti->aseq, proti->ncont, proti->len,
		  protj->Cont_map, protj->aseq, protj->ncont, protj->len);

    Write_identity(it1, proti, protj, ali_ij,
		   id_aa[it1], id_sup[it1], shift[it1],
		   neigh_ali[it1], neigh_noali[it1], neigh_noali_aaid[it1],
		   ali_cont_ct[it1], id_cont_ct[it1],
		   ali_cont_sup[it1], id_cont_sup[it1],
		   ali_cont_aaid[it1], id_cont_aaid[it1]);
  }
}
void Sum_scores(double *nali_sum, double *SI_sum, double *TM_sum,
		double *CO_sum, double *PC_sum,
		float *nali, float *SI, float *TM,
		float *CO, float *PC, int it){
  nali_sum[it]+=nali[it];
  SI_sum[it]+=SI[it];
  TM_sum[it]+=TM[it];
  CO_sum[it]+=CO[it];
  PC_sum[it]+=PC[it];
}


void Print_scores(char *out, double *s1, char *what){
  char line[800], tmp[30];
  sprintf(line, "Mean %15s %s: %.4f", what, ali_name[0], s1[0]/npair_pdb);
  if(ALI_SS){
    sprintf(tmp," sec.str_ali: %.4f", s1[code[1]]/npair_pdb); strcat(line, tmp);
  }
  if(ALI_TM){
    sprintf(tmp," TM_ali: %.4f", s1[code[2]]/npair_pdb); strcat(line, tmp);
  }
  if(ALI_CO){
    sprintf(tmp," CO_ali: %.4f", s1[code[3]]/npair_pdb); strcat(line, tmp);
  }
  if(ALI_PC){
    sprintf(tmp," PC_ali: %.4f", s1[code[4]]/npair_pdb); strcat(line, tmp);
    int it=PTYPE;
    for(int i=0; i<NMSA_act; i++){
      sprintf(tmp," %s: %.4f",ali_name[it],s1[it]/npair_seq); it++;
      strcat(line, tmp);
    }
  }
  strcat(line, "\n");
  printf("%s", line);
  strcat(out, line);
}

// Fortran routines, you may need to change the underscores
// depending on fortran compiler
#define SIMAT simat_c__
#define MAXSUP maxsup_c__


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include "allocate.h"

// --- protein descriptions

#define NAMAX 20 // Maximum number of atoms in a residue
#define MAXRES 5000
#define MAXATOM 30
#define CHARPDB 10
#define NCHAR 200

struct atom{
  float r[3];
  char name[4];
};
struct residue{
  int n_atom;
  struct atom atom[MAXATOM];
  char aa;
  char label[6];
};
struct sec_el{
  int ini, len;
  char type;
};

struct protein{
  int len;
  int nca;
  char name[CHARPDB];
  char name_file[100];
  char chain[80];
  char *ss;
  char *aseq;
  int *nseq;
  int *n_atom;
  char **pdbres;
  float **xca;
  float **vec;
  struct atom **res_atom;
  //
  char exp_meth;
  int N_cont, n_sec_el;
  char  *seq2;
  short *ss_num;
  short **Cont_map;
  int *ncont;
  float *EC;
  struct sec_el *sec_el;
};

// Amino acid names
#define AANAME1 "ACDEFGHIKLMNPQRSTVWYX"
#define AANAME3 "ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYRXXX"


char *modres[20];
int  n_modres[20];

// 1mammoth.c
void mammoth(struct protein *prot1, int *iv1, struct protein *prot2, int *iv2,
	     char *match1, char *match3, int *nali,
	     int use_long,  int *norm,float *zscore,
	     float *lne, double *evalue, float *scorel,
	     int *nsup_1, float *psi_1, float *rms_1, 
	     int *nsup_3, float *psi_3, float *rms_3,
	     float *rotation, float *x1_ave, float *x2_ave,
	     float *tm_score, float *seq_sim);

// align.c
int align(int report, float **smap, float gapi, float gape, int maxl,
	  int nres_1, int nres_2, int *iv1, int *iv2, char *match,
	  float SIMTHR);

//align3d.c
int Simat_3D(float **smat3d, float **smat1d, float *rotation,
	     float **xca1, float *x1_ave, int nres_1,
	     float **xca2, float *x2_ave, int nres_2,
	     int *iv1, int *iv2, char *match3, int nali);

// aux_mammoth.c
int Check_file_in(char *filename, int exit);
void Modres();
char **getpdb_names(char *list_pdb, int *n_pdb, char **chain_lab);

// output.c
int Output(FILE *file_out, struct protein *prot1, struct protein *prot2,
	   int terse, int generate_pdb, float z_thresh,
	   int *iv1, int *iv2, int nvnew,
	   char *match1, int nsup_1, float psi_1, float rms_1,
	   char *match3, int nsup_3, float psi_3, float rms_3,
	   float zscore, float lne, double evalue, float scorel, int norm,
	   float tm_score, float seq_sim);

// read_pdb_mammoth.c
int Read_pdb(struct protein *prot, char *filename, char *chain_to_read);

// Read_command_line_mammoth.c

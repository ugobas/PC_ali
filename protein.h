struct atom{
  float r[3];
  char name[4];
};
struct residue{
  int n_atom;
  struct atom *atom;
  //struct atom atom[MAXATOM]
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
  int natoms;
  char name_file[20];
  char code[30];
  char domname[40];
  char domain[40];
  char chain;
  //char *ss;
  char *aseq;
  short *seq_bs;  // blosum code
  short *seq1;    // Alphabetic one-letter code
  short *ss3;     // sec.str. in three letters code
  short *n_atom;  // Number of atoms per residue
  float *max_sim; // Sim. score of id. prot. 0=ali 1=blosum 2=ss 3=TM 4=cont
  char **pdbres;
  short *seqres; // Numbering according to seqres
  char *seqr;     // One-letter code for seqres residues
  float *xca_rot; // alpha carbon coordinates in optimal superimposition
  float **vec;
  struct atom **res_atom;
  char exp_meth;
  //
  int N_cont, n_sec_el;
  //char  *seq2;
  short *ss_num;
  short **Cont_map;
  int *ncont;
  float *EC;
  struct sec_el *sec_el;
};

// Amino acid names
#define AANAME1 "ACDEFGHIKLMNPQRSTVWYX"
#define AANAME3 "ALACYSASPGLUPHEGLYHISILELYSLEUMETASNPROGLNARGSERTHRVALTRPTYRXXX"

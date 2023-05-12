float TM_score(float **d2, float *d02, float **rot_out, float *shift_out, 
	       int *ali, float ltar,
	       struct protein *prot1, struct protein *prot2, int verbose);
float TM_score_mult(float **TM_all, int **msa, int N_ali, int nprot,
		    struct protein **prot_p, int verbose);

void Align_TM(int *ali_new, float *d2min1, float **d2, float d02, int *ali,
	      int n1, int n2);
float TM_fast(float **d2_out, float d02, int *ali, int ltar, float *d2min1,
	      struct protein *prot1, struct protein *prot2);
void Examine_neighbors(float **d2, int *ali, float d02,
		       int *shift, int *id_3D,
		       int *neigh_ali, int *neigh_noali, int *neigh_noali_aaid,
		       char *seq1, int n1, char *seq2, int n2);
void Test_Rot(float **d2, int i, int j, float **Rot, float *Shift,
	      float *xca1_store, int len1, float *xca2_store, int len2);

void Copy_coord(float *xca, float *xca_store, int len);
void Rot_shift(float *xca_i, int len, float **rot, float *shift);

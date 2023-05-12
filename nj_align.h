int NJ_align(int **msa, int *N_ali_max,
	     int ***Ali_pair, float ****Rot_pair, float ***Shift_pair,
	     struct protein **prot_p, char *name_out, 
	     float **div_in, int n, int AVE_LINK);
int Align_PC_NW_mult(int **msa, int *len_msa, int *N_ali_max, int N_seq,
		     struct protein **prot_p,
		     int ***Ali_pair,float ****Rot_pair,float ***Shift_pair,
		     int *index1, int n1, int *index2, int n2);
// Align_pair[i][j][l]= residue in protein j aligned to res l of protein i
int Align_PC_NW_pair(int *ali_PC, float **d2, int **nc,
		     struct protein *proti, struct protein *protj);
void Set_scores(struct protein *prot);

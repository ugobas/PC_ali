void Clique_MSA(int **msa, int *L_msa,
		int N_ali_max,
		int ***Ali_pair, // Ali_pair[i][j][site_i]=site_j i<j
		int *Lprot,      // Length of protein i
		int Nprot,       // Number of protein
		char **Seq,      // All sequences
		char **name_seq, // Names of sequences
		char *name_in,
		int L_msa_ini);   // Length of the input MSA
void Print_MSA(int **msa_new, char *name_msa, int Nprot, int L_msa,
	       int *Lprot, char **Seq, char **name_seq);

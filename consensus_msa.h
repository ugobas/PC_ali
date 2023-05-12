void Consensus_MSA(int ***Ali_pair, // Ali_pair[i][j][site_i]=site_j i<j
		   int **PC_opt_pair,
		   int *Lprot,      // Length of protein i
		   int Nprot,       // Number of proteins
		   int **MSA,       // MSA[i][column k]=site_i
		   int L_msa,       // Length of the MSA
		   char **Seq,      // All sequences
		   char **name_seq, // Names of sequences
		   char *name_in);

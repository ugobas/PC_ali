int Read_structures(struct protein **prots1, char *PDB_PATH, char *FILE_LIST,
		    char *EXT_CM, char *EXT_SEQ, int CONT_DEF,
		    char CONT_TYPE, float CONT_THR, int IJ_MIN);
int Read_PDB_compress(struct protein **prot,
		      char *pdbid, char *chain,
		      char *PDB_PATH, char *PDB_EXT);
int Read_pdb(char *filename, struct protein **prot, char *chain_to_read);
int Select_domain(struct protein *prot, 
		  int *ini_frag, int *end_frag, int nfrag);

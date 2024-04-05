
int Align_seq(int *Prot_ali, int N_ali, char *Prot_seq, char *PDB_seq, int n2);
int Pair_ali(int *ali, int N_ali, int *ali1, int *ali2);
float Seq_identity(char *seq_i, char *seq_j, int N_ali);
float Contact_overlap(int *ali,
		      short **Cont_mat_1, int N1,
		      short **Cont_mat_2, int N2);
float Test_contacts(int *ali, int *id_3D, float c_ave,
		    int *ali_cont_ct,  int *id_cont_ct,
		    int *ali_cont_sup, int *id_cont_sup,
		    int *ali_cont_aaid,int *id_cont_aaid, 
		    short **Cont_mat_1, char *seq1, int *ncont1, int N1,
		    short **Cont_mat_2, char *seq2, int *ncont2, int N2);
void Shared_contacts(int **nc, int *ali,
		     short **Cont_mat_1, int n1,
		     short **Cont_mat_2, int n2);
void Shared_contacts_noali(int **nc,
			   short **Cont_mat_1, int n1,
			   short **Cont_mat_2, int n2);
void Align_CO(int *ali_CO, int *ali, int **nc,
	      short **Cont_map1, int n1,
	      short **Cont_map2, int n2);
float Seqid(int *diff, int *ali_12, int *ali_aa,
	    char *seq1, int N1, char *seq2, int N2);
char Get_compression(char *file_name);
int Find_name(char *name, char **names, int N, int Np);
int Change_ext(char *name_out, char *name_in, char *ext);
void Align_neighbors(int *ali_new, int *ali, float thr,
		     int *neighb1, float *score1, int n1,
		     int *neighb2, float *score2, int n2);
void Invert_ali(int *ali2, int n2, int *ali, int n1);

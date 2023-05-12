struct treenode{
  struct treenode *offs1, *offs2;
  struct treenode *parent;
};

struct treenode *Neighbor_Joining(float **div, int n);
int ***Build_outgroups(int ***N_out, struct treenode *node, int n);
int ***Outgroups_NJ(int ***N_out, float **div, int n);
int Single_linkage(int ***elements, int **n_ele, int **clus,
		   float **diff, int n, float thr);

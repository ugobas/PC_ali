#include "protein.h"
#include "tm_score.h"
#include "McLachlan.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define AA3 "ALA GLU GLN ASP ASN LEU GLY LYS SER VAL ARG THR PRO ILE MET PHE TYR CYS TRP HIS "
#define AA1 "AEQDNLGKSVRTPIMFYCWH"

void Name3(char *aaname3, char aa1);
extern int Copy_CA(struct atom *atom1, float *xca, int n_atom);
void Print_pdb(char *name_in, struct protein **prot_p, int n)
{
  char nameout[100]; sprintf(nameout, "%s.pdb", name_in);
  FILE *file_out=fopen(nameout, "w");
  printf("Printing %d superimposed structures in %s\n", n, nameout);
  fprintf(file_out, "REMARK %d structures superimposed after multiple "
	  "alignment by maximizing the multiple TM score\n", n);

  char aaname3[10]; int i; for(i=0; i<4; i++)aaname3[i]='\0';
  float shift[3], *rot[3];
  for(i=0; i<3; i++)rot[i]=malloc(3*sizeof(float)); 
  for(int ip=0; ip<n; ip++){
    struct protein *prot=prot_p[ip];
    char chain=prot->chain;
    int L=prot->len;
    fprintf(file_out, "MODEL %d %s L=%d\n", i+1, prot->name, L);
    float xca_pdb[3*L], *xca=xca_pdb, w[L]; int nca=0;
    for(int i=0; i<L; i++){
      if(Copy_CA(prot->res_atom[i], xca, prot->n_atom[i])==0){
	printf("WARNING, %d atoms found at res %d%c of %d\n",
	       prot->n_atom[i], i, prot->aseq[i], L); continue;
      }
      xca+=3; w[i]=1; nca++;
    }
    // Align to rotated structures
    rmsd_mclachlan_f(rot, shift, prot->xca_rot, xca_pdb, w, nca);
    // print 
    int na=0;
    for(int i=0; i<L; i++){
      Name3(aaname3,prot->aseq[i]);
      struct atom *atomi=prot->res_atom[i];
      for(int j=0; j<prot->n_atom[i]; j++){
	Rot_shift(atomi->r, 1, rot, shift);
	if(na < 99999)na++;
	fprintf(file_out,
		"ATOM  %5d  %3s %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
		na, atomi->name, aaname3, chain, prot->pdbres[i],
		atomi->r[0],atomi->r[1],atomi->r[2],1.0, 80.0);
	atomi++;
      }
    }
    fprintf(file_out, "TER\n");
    fprintf(file_out, "ENDMDL\n");
  }
  fclose(file_out);
  for(i=0; i<3; i++)free(rot[i]);
}

void Name3(char *aaname3, char aa1){
  int a; for(a=0; a<20; a++){if(aa1==AA1[a])break;}
  if(a==20){printf("WARNING, aa %c not found\n", aa1); a=0;}
  a*=4; for(int j=0; j<3; j++){aaname3[j]=AA3[a]; a++;}
  aaname3[3]='\0';
}

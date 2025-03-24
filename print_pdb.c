#include "protein.h"
#include "tm_score.h"
#include "McLachlan.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
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
    fprintf(file_out, "MODEL %d %s L=%d\n", i+1, prot->domname, L);
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

float Compute_RMSD(struct protein **prot_p, int ci, int cj,
		   int **msa, int L_msa)
{
  float *xca_i=prot_p[ci]->xca_rot;
  float *xca_j=prot_p[cj]->xca_rot;
  int Li=prot_p[ci]->len;
  int Lj=prot_p[cj]->len;
  int lmin=Li; if(Lj<lmin)lmin=Lj;

  float TM_coeff=1.24, TM_offset=1.8, L_offset=15;
  float d02=TM_coeff*pow(lmin-L_offset, 1./3)-TM_offset; //lmin
  //if(d02<=0)return(0);
  d02*=d02;
  int n=0; float rmsd=0;
  for(int k=0; k<L_msa; k++){
    int ki=msa[ci][k], kj=msa[cj][k];
    if(ki<0 || kj<0)continue;
    float *x1=xca_i+3*ki, *x2=xca_j+3*kj;
    float dx=x1[0]-x2[0], dy=x1[1]-x2[1], dz=x1[2]-x2[2];
    float dd=dx*dx+dy*dy+dz*dz;
    //if(dd>d02)continue; // Compute only for well superimposed
    if(dd>d02){dd=d02;}; // Min between dd and d02
    rmsd+=dd;
    n++;
  }
  if(n)rmsd/=n;
  return(sqrt(rmsd));
}

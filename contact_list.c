#define CONT_MAX 30000
int ini_cont=0;
float C_THR, CONT_THR_2;

#include "protein.h"
#include "cont_list.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/******** Internal routines *********************/
struct atom *Find_atom(struct atom *atom, char *ATOM, int n);
int Contact(struct atom *atom_i, int ni, struct atom *atom_j, int nj,
	    char *ATOM);
int Contact_all(struct atom *atom_i, int ni, struct atom *atom_j, int nj);


/*************************** Main routine *********************/
int Compute_contact_list(struct protein *prot, char cont_type,
			 float cont_thr, int IJ_MIN)
{
  // INPUT
  int N_res=prot->len;
  struct atom **atoms=prot->res_atom, *atom_i;
  short *n_atom=prot->n_atom, ni;

  // OUTPUT
  int NC=0;
  short **Cont_list=malloc(N_res*sizeof(short *));
  int *n_cont=malloc(N_res*sizeof(int));

  // Auxiliary
  int i_res, j_res, i, contact=0;
  short res1[CONT_MAX], res2[CONT_MAX];
  int *num_cont=malloc(N_res*sizeof(int));

  if(ini_cont==0){
    C_THR=cont_thr;
    CONT_THR_2=cont_thr*cont_thr;
    printf("Contact type: %c Threshold: %.2f A |i-j|>%d\n",
	   cont_type,cont_thr,IJ_MIN);
    ini_cont=1;
  }

  for(i_res=0; i_res<N_res; i_res++){
    num_cont[i_res]=0;
    atom_i=atoms[i_res]; ni=n_atom[i_res];
    for(j_res=i_res+IJ_MIN; j_res< N_res; j_res++){
      if(cont_type=='c'){
	contact=Contact_all(atom_i, ni, atoms[j_res], n_atom[j_res]);
      }else if(cont_type=='b'){
	contact=Contact(atom_i, ni, atoms[j_res], n_atom[j_res], "CB");
      }else if(cont_type=='a'){
	contact=Contact(atom_i, ni, atoms[j_res], n_atom[j_res], "CA");
      }
      if(contact){
	if(NC>=CONT_MAX){
	  printf("ERROR, more than %d contacts found\n", NC); exit(8);
	}
	res1[NC]=i_res; res2[NC]=j_res; NC++;
	num_cont[i_res]++; 
      }
    }
  }

  for(i_res=0; i_res<N_res; i_res++){
    Cont_list[i_res]=malloc((num_cont[i_res]+1)*sizeof(short));
    Cont_list[i_res][num_cont[i_res]]=-1;
    num_cont[i_res]=0;
    n_cont[i_res]=0;
  }
  for(i=0; i<NC; i++){
    i_res=res1[i]; j_res=res2[i];
    n_cont[i_res]++; n_cont[j_res]++;
    Cont_list[i_res][num_cont[i_res]]=j_res;
    num_cont[i_res]++;
  }

  prot->N_cont=NC;
  prot->ncont=n_cont;
  prot->Cont_map=Cont_list;
  free(num_cont);
  return(NC);
}


int Contact_all(struct atom *atom_i, int ni, struct atom *atom_j, int nj)
{
  struct atom *atom1=atom_i, *atom2;
  float dx, dy, dz;
  int i, j;

  for(i=0; i<ni; i++){
    atom2=atom_j;
    for(j=0; j<nj; j++){
      dx=(atom1->r[0]-atom2->r[0]); if(fabs(dx)>C_THR) goto new;
      dy=(atom1->r[1]-atom2->r[1]); if(fabs(dy)>C_THR) goto new;
      dz=(atom1->r[2]-atom2->r[2]); if(fabs(dz)>C_THR) goto new;
      if((dx*dx+dy*dy+dz*dz)<=CONT_THR_2) return(1);
    new:
      atom2++;
    }
    atom1++;
  }
  return(0);
}

int Contact(struct atom *atom_i, int ni, struct atom *atom_j, int nj,
	    char *ATOM)
{
  struct atom *atom1=Find_atom(atom_i, ATOM, ni);
  struct atom *atom2=Find_atom(atom_j, ATOM, nj);
  float dx, dy, dz;

  if((atom1==NULL)||(atom2==NULL))return(0);
  dx=(atom1->r[0]-atom2->r[0]); if(fabs(dx)>C_THR) return(0);
  dy=(atom1->r[1]-atom2->r[1]); if(fabs(dy)>C_THR) return(0);
  dz=(atom1->r[2]-atom2->r[2]); if(fabs(dz)>C_THR) return(0);
  if((dx*dx+dy*dy+dz*dz)<=CONT_THR_2)return(1);
  return(0);
}

struct atom *Find_atom(struct atom *atom, char *ATOM, int n)
{
  struct atom *atom1=atom; int i;
  for(i=0; i<n; i++){
    if(strncmp(atom1->name, ATOM, 2)==0)return(atom1);
    atom1++;
  }
  // Not found, return CA
  atom1=atom;
  for(i=0; i<n; i++){
    if(strncmp(atom1->name, "CA", 2)==0)return(atom1);
    atom1++;
  }
  return(NULL);
}

#include "protein.h"
#include "allocate.h"
#include "read_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define SEC_STR_MAX 500      // Max number of sec. str.elements
#define MAXRES 10000
#define MAXATOM 100
#define CHARPDB 20

// Secondary structure types: " HE"
#include "blosum62.h"
#include "secstrsim.h"
//
//char AA_Blosum[22]= "ARNDCQEGHILKMFPSTWYVX*";
//char SS_code[]=" HE-";

int Code_AA(char res);
int Code_AA_2(char aseq, char *aacode, int n);
char Get_aaname(char *aaname);
void Set_prot_name(char *name, char *filename, char chain);
int Write_protein(struct protein *prot, struct residue *res, int nres,
		  char *filename, char chain_to_read, char exp_meth);
int Copy_CA(struct atom *atom1, float *xca, int n_atom);
void Set_CA_vectors(struct protein *prot);

struct secondary{
  char ini_res[6], end_res[6];
  char type, chain;
};
static int Read_sec_str(struct secondary *sec_ele, int *N_sec_ele,
		 char type, char *string);
static void Copy_word(char *word, char *origin, int ini, int end);
static int Secondary_structure(short *ss3, int N_res, char chain, char **pdbres,
			       struct secondary *sec_ele, int N_sec_ele);
int Ini_Modres=0;
void Modres();
char *modres[20];
int  n_modres[20];

/***********************************************************/

int Read_pdb(char *filename, struct protein **prot, char *chain_to_read)
{

  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist, exiting\n", filename);
    return(0);
  }
  printf("Reading PDB: %s\n",filename);
  if((*chain_to_read=='_')||(*chain_to_read==' '))*chain_to_read='\0';

  char string[300];
  int nres=0, n_atom=0;
  int alternative=0, alt_check=0, het=0, exo=1, n_exo=0;
  char alt_1=' ', alt, icode_old=' ',  icode, chain, exp_meth=' ';
  char res_label[6], res_old[6]="xxxx\0";
  float x, y, z;

  struct atom *atom;
  struct residue res[MAXRES], *res_i=res;
  for(int i=0; i<MAXRES; i++)res[i].atom=NULL;
  struct secondary sec_ele[SEC_STR_MAX];
  int N_sec_ele=0;

  while(fgets(string, sizeof(string), file_in)!=NULL){

    if(res_i->atom==NULL)res_i->atom=malloc(MAXATOM*sizeof(struct atom));
    if(strncmp(string, "ATOM", 4)==0){
      // Standard residue or DNA basis
      het=0;
    }else if(strncmp(string, "HETATM", 6)==0){
      // Cofactor or exotic residue
      het=1;
    }else if((strncmp(string,"TER",3)==0)&&(nres>0)&&(*chain_to_read!='*')){
      break;
    }else if(strncmp(string, "ENDMDL", 6)==0){
      printf("WARNING: pdb entry contains multiple models\n");
      printf("WARNING: only the first model will be used\n");
      break;
    }else if(strncmp(string,"EXPDTA", 6)==0){
      if(strncmp(string+10, "NMR", 3)==0){exp_meth='N';}
      else{exp_meth='X';}
      continue;
    }else if(strncmp(string,"HELIX ", 6)==0){
      Read_sec_str(sec_ele, &N_sec_ele, 'H', string); continue;
    }else if(strncmp(string,"SHEET ", 6)==0){
      Read_sec_str(sec_ele, &N_sec_ele, 'E', string); continue;
    }else{
      continue; 
    }
 
    // Check chain
    chain=string[21];
    if((nres<=0)&&((*chain_to_read)=='\0'))*chain_to_read=chain; 
    if((*chain_to_read)=='.'){chain='.';}
    if(((*chain_to_read)!='*')&&((*chain_to_read)!=chain)){
      if(nres>0){break;}else{continue;}
    }
 
    // Omit hydrogen atoms
    if(string[13]=='H')continue;

    if(het){
      //   if(HETATM), check if amino acid (N-CA-...)
      if(n_exo==0){
	if(strncmp(string+13, "N ", 2)!=0){exo=0; continue;}
	else{n_exo=1; exo=1;}
      }else if(exo==0){
	continue;
      }else if(n_exo==1){
	if(strncmp(string+13, "CA", 2)!=0){exo=0; n_atom=0; continue;}
	else{n_exo=2;}
      }
    }

    // Alternative configurations eliminated (indicated as ALT2)
    if((strncmp(string+72, "ALT", 3)==0)&&(string[75]!='1'))continue;
     
    // Alternative configurations eliminated (indicated as ARES)
    alt=string[16];
    if(alt!=' '){
      if(alt_1==' ')alt_1=alt;
      if(alt!=alt_1)continue;
    }

    // Read
    icode=string[26]; //string[26]=' ';
    sscanf(string+22, "%s%f%f%f", res_label, &x, &y, &z);

    // Alternative configurations eliminated (indicated as insertion code)
    if((icode!=icode_old)&&(strncmp(res_label,res_old,4)==0)){
      if(alternative)continue;
	if(alt_check==0){
	  struct atom atom_old=res_i->atom[0];
	  float dx=x-atom_old.r[0], dy=y-atom_old.r[1], dz=z-atom_old.r[2];
	  alt_check=1;
	  if((dx*dx+dy*dy+dz*dz)<1){
	    alternative=1; continue;
	  }
	}
    }

    // New residue ?
  if(strncmp(res_label,res_old,5)!=0){
      // Write old residue
      if(nres>=0)res_i->n_atom=n_atom;
      // Write new residue
      nres++; res_i++;
      if(res_i->atom==NULL)
	res_i->atom=malloc(MAXATOM*sizeof(struct atom));
      n_atom=0;
      res_i->aa=Get_aaname(string+17);
      strcpy(res_i->label, res_label);
      // Store
      strcpy(res_old, res_label);
    }

    if(nres >= MAXRES){
      printf("ERROR, too many residues (more than %d)\n", MAXRES);
      exit(8);
    }

    // Write atom
    atom=res_i->atom+n_atom;
    atom->r[0]=x; atom->r[1]=y; atom->r[2]=z;
    atom->name[0]=string[13]; atom->name[1]=string[14];
    atom->name[2]=string[15]; atom->name[3]='\0';
    n_atom++;
    
    if(n_atom >= MAXATOM){
      printf("ERROR, too many atoms in residue %d %c (>%d)\n",
    	     nres, res_i->aa, n_atom); exit(8);
    }

  }
  fclose(file_in);

  if((n_atom>0)&&((het==0)||(n_exo>=2))){res_i->n_atom=n_atom; nres++;}

  if(nres<=0){
    printf("ERROR, no amino acid found in pdb file %s chain %c\n",
	   filename, *chain_to_read); return(0);
  }
  printf("%d amino acid found in pdb %s chain %c\n",
	 nres, filename, *chain_to_read);

  // Write protein
  Write_protein(*prot, res, nres, filename, *chain_to_read, exp_meth);
  Secondary_structure((*prot)->ss3, (*prot)->len, *chain_to_read, //nres
		      (*prot)->pdbres, sec_ele, N_sec_ele);
  //printf("%d secondary structure elements found\n", N_sec_ele);

  for(int i=0; i<nres; i++)if(res[i].atom)free(res[i].atom);

  return(nres);
}

int Write_protein(struct protein *prot, struct residue *res, int nres,
		  char *filename, char chain_to_read, char exp_meth)
{
  struct atom *atom1, *atom2; 
  int i, j, nca=0, natoms=0;
  for(i=0; i<nres; i++)natoms+=res[i].n_atom;
  prot->natoms=natoms;

  Set_prot_name(prot->name, filename, chain_to_read);
  prot->exp_meth=exp_meth;
  //prot->ss=malloc(nres*sizeof(char));
  prot->ss3= malloc(nres*sizeof(short));
  prot->aseq=malloc(nres*sizeof(char));
  prot->seq1=NULL;
  prot->seq_bs=malloc(nres*sizeof(short));
  prot->n_atom=malloc(nres*sizeof(short));
  prot->xca_rot=malloc((3*nres)*sizeof(float));
  prot->vec=NULL;   //=Allocate_mat2_f_fortran(nres,3);
  prot->res_atom=malloc(nres*sizeof(struct atom *));
  prot->pdbres=malloc(nres*sizeof(char *));
  prot->max_sim=NULL;

  //prot->xall=malloc((3*atoms)*sizeof(float));
  //float *xall=prot->xall;

  float *xca_rot=prot->xca_rot;
  for(i=0; i<nres; i++){
    int n_atom=res[i].n_atom;

    // Copy residue only if CA is found:
    if(Copy_CA(res[i].atom, xca_rot, n_atom)==0){
      printf("WARNING, %d atoms found at res %d%c of %d\n",
	     n_atom, i, res[i].aa, nres); continue;
    }

    prot->aseq[nca]=res[i].aa;
    prot->n_atom[nca]=res[i].n_atom;
    prot->pdbres[nca]=malloc(6*sizeof(char));
    strcpy(prot->pdbres[nca], res[i].label);
    if(prot->seq1)prot->seq1[nca]=Code_AA_2(res[i].aa, AANAME1, 20);
    if(prot->seq_bs)prot->seq_bs[nca]=Code_AA_2(res[i].aa, AA_Blosum, 22);

    // Write atoms
    prot->res_atom[nca]=malloc(n_atom*sizeof(struct atom));
    atom1=prot->res_atom[nca]; atom2=res[i].atom;
    for(j=0; j<n_atom; j++){
      *atom1=*atom2; atom1++; atom2++;
    }
    nca++; xca_rot+=3;
  }
  prot->nca=nca;
  prot->len=nca;
  if(nca==0){
    printf("ERROR, no CA atoms found in pdb file %s\n", filename); exit(8);
  }
  if(prot->vec)Set_CA_vectors(prot);
  //if(ASSIGN_SS)Assign_ss(prots+i_prot)

  return(0);
}

int Copy_CA(struct atom *atom_first, float *xca, int n_atom)
{
  if(n_atom<=0){
    printf("WARNING in Copy_CA, %d atoms", n_atom);
    if(n_atom)printf(" first= %s", atom_first->name);
    printf("\n");
    return(0);
  }
  struct atom *atom=atom_first;
  for(int i=0; i<n_atom; i++){
    if((atom->name[0]=='C')&&(atom->name[1]=='A'))goto found;
    atom++;
  }
  printf("WARNING, CA not found using first atom\n");
  atom=atom_first;
  
 found: // CA found
  xca[0]=atom->r[0];
  xca[1]=atom->r[1];
  xca[2]=atom->r[2];
  if(isnan(xca[0]) || isnan(xca[1]) || isnan(xca[2])){
    printf("ERROR in atom %s r= %.2g %.2g %.2g n_atom= %d\n",
	   atom->name, atom->r[0], atom->r[1], atom->r[2], n_atom);
    exit(8);
  }

  return(1);
}

char Get_aaname(char *aaname){
  int i, j;
  if(Ini_Modres==0){Modres(); Ini_Modres=1;}
  for(i=0; i<20; i++){
    if(strncmp(aaname, AANAME3+3*i, 3)==0)return(AANAME1[i]);
  }
  for(i=0; i<20; i++){
    for(j=0; j<n_modres[i]; j++){
      if(strncmp(aaname, modres[i]+3*j, 3)==0)return(AANAME1[i]);
    }
  }
  return('X');
}

void Set_prot_name(char *name, char *filename, char chain){
  int i=0; char *ptr=filename, tmp[400];
  //if(name[0]!='\0')return;
  // copy name without directory path
  if((*ptr=='.')&&(*(ptr+1)=='/'))ptr+=2;
  while((*ptr!='\n')&&(*ptr!='\0')&&(*ptr!='.')){
    if(*ptr=='/'){ptr++; i=0;}
    tmp[i]=*ptr; ptr++; i++;
  }
  // Add chain label
  if(i<CHARPDB){if(chain!=' '){tmp[i]=chain;}else{tmp[i]='_';} i++;}
  tmp[i]='\0';
  //printf("%d characters in name\n", i);
  if(i>=CHARPDB){printf("WARNING, too many characters in prot name\n");}
  else{strcpy(name, tmp);}
}

/* Numeric codes */
int Code_AA(char res){
  int i; 
  for(i=0; i<20; i++)if(res==AANAME1[i])return(i);
  if((res!='-')&&(res!='.')&&(res!='*'))
    printf("Warning, wrong aa type %c\n", res);
  return(20);
}


void Set_CA_vectors(struct protein *prot){
  int nv=prot->len-1, i, k;
  float **vec=prot->vec, *xca=prot->xca_rot, *xca1=xca+3;

  for(i=0; i<nv; i++){
    float vmod=0;
    float *v=vec[i];
    for(k=0; k<3; k++){v[k]=xca1[k]-xca[k]; vmod+=v[k]*v[k];}
    vmod=sqrt(vmod);
    for(k=0; k<3; k++)v[k]/=vmod;
    xca=xca1; xca1+=3;
  }
}

int Read_sec_str(struct secondary *sec_ele, int *N_sec_ele,
		 char type, char *string)
{
  if((*N_sec_ele)>=SEC_STR_MAX){
    printf("WARNING Too many secondary structure elements (more than %d)\n",
	   *N_sec_ele); return(-1);
  }
  sec_ele[*N_sec_ele].type=type;

  /* initial and final residues */
  if(type=='H'){
    sec_ele[*N_sec_ele].chain=string[19];
    Copy_word(sec_ele[*N_sec_ele].ini_res, string, 21, 25);
  }else if(type=='E'){
    sec_ele[*N_sec_ele].chain=string[21];
    Copy_word(sec_ele[*N_sec_ele].ini_res, string, 22, 26);
  }
  Copy_word(sec_ele[*N_sec_ele].end_res, string, 33, 37); 
  (*N_sec_ele)++;
  return(0);
}

void Copy_word(char *word, char *origin, int ini, int end){
  char *w=word, *o=origin+ini; int i;
  for(i=ini; i<=end; i++){*w='\0'; w++;}
  w=word;
  for(i=ini; i<=end; i++){
    if((*o!=' ')&&(*o!='\0')){*w=*o; w++;} o++;
  }
}

int Secondary_structure(short *ss3, int N_res, char chain, char **pdbres,
			struct secondary *sec_ele, int N_sec_ele)
{
  int i_sec, i_res;
 
  printf("Storing %d secondary structure elements\n",N_sec_ele);
  for(i_res=0; i_res<N_res; i_res++)ss3[i_res]=0;
  for(i_sec=0; i_sec<N_sec_ele; i_sec++){
    if(sec_ele[i_sec].chain!=chain)continue;
    for(i_res=0; i_res<N_res; i_res++){
      if(strncmp(pdbres[i_res], sec_ele[i_sec].ini_res, 5)==0)break;
    }
    //printf("%d %c %s %s  %s\n", i_sec,sec_ele[i_sec].type,
    //	   sec_ele[i_sec].ini_res,sec_ele[i_sec].end_res, pdbres[i_res]);
    while(i_res<N_res){
      if(sec_ele[i_sec].type=='H'){ss3[i_res]=1;}
      else if(sec_ele[i_sec].type=='E'){ss3[i_res]=2;}
      if(strncmp(pdbres[i_res], sec_ele[i_sec].end_res, 5)==0)break;
      i_res++;
    }
  }
  for(i_res=0; i_res<N_res; i_res++)printf("%c",SS_code[ss3[i_res]]);
  printf("\n");
  //exit(8);
  return(N_sec_ele);
}

void Modres(){
  int i;
  for(i=0; i<20; i++)modres[i]=malloc(100*sizeof(char));
  strcpy(modres[0], "CSDMDOORNDBZNALLALHACAYANCB"); // Ala
  strcpy(modres[1], // Cys 
	 "CYGCY3CSUCSPCYMCASCSBCSRCMECMHCSSCSXCSWCSOALSSMCCEAOCSCCSSCYYCM");
  strcpy(modres[2], "BFDSNNASABHDSUI"); // Asp
  strcpy(modres[3], "PCAAR4GLQGMA");    // Glu 
  strcpy(modres[4], "DPN"); // Phe
  strcpy(modres[5], "GL3CHGGLZACYFGL"); //  Gly 
  strcpy(modres[6], "MHSNEPHIC"); //  His
  strcpy(modres[7], "IIL"); // Ile
  strcpy(modres[8], "LLPKCXLYZALYLCXMCL"); // Lys 
  strcpy(modres[9], "LEFDLE"); //  Leu 
  strcpy(modres[10], "MSEFMEMMEMSOFORNRQCH6"); // Met
  strcpy(modres[11], "IASMEN"); // Asn
  strcpy(modres[12], "HYPDPR"); // Pro 
  strcpy(modres[13], "MGN5HP"); // Gln 
  strcpy(modres[14], "AGMBORARMACLDAR"); // Arg
  strcpy(modres[15], "OSEDSNHSLDHLSACSEP"); // Ser
  strcpy(modres[16], "TPO"); //  Thr
  strcpy(modres[17], "DVA"); // Val
  strcpy(modres[18], "TROTRQTRWTRN"); //  Trp
  strcpy(modres[19], "TPQPTRSTYYOFTYITYCTYS"); // Tyr
  for(i=0; i<20; i++){
    int k=0; while(modres[i][k]!='\0')k++; n_modres[i]=k/3;
    //printf("%c %d\n", AANAME1[i], n_modres[i]);
  }
}  

int Code_AA_2(char aseq, char *aacode, int n){
  for(int i=0; i<n; i++)if(aseq==aacode[i])return(i);
  if(aseq=='X')return(0);
  printf("WARNING, aa %c not known\n", aseq);
  return(20);
}


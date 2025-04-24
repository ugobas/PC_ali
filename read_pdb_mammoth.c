#include "protein.h"
#include "allocate.h"
#include "read_structures.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#define SEC_STR_MAX 5000      // Max number of sec. str.elements
#define MAXRES 10000
#define MAXATOM 100
#define CHARPDB 20

// Secondary structure types: " HE"
#include "blosum62.h"
#include "secstrsim.h"
//
//char AA_Blosum[22]= "ARNDCQEGHILKMFPSTWYVX*";
//char SS_code[]=" HE-";

// Modified residues
#define NEXO 400
  int n_exo=0;
char *res_exo[NEXO], *res_std[NEXO];

int Read_modres(char **res_exo, char **res_std, char *string, int *n_exo);
char Het_res(char *res_type_old, char **res_exo, char **res_std, int n_exo);
int Code_AA(char res);
int Code_AA_2(char aseq, char *aacode, int n);
int Get_number(char *pdbres);
char Get_aaname(char *aaname);
char Code_3_1(char *res);
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
static short *Read_seqres(struct protein *prot, int nres,
			  char *filename, char chain_to_read);
static void Get_seqres(char *string, char *seqres, int *n_seqres);

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
    }else if(strncmp(string,"MODRES", 6)==0){
      Read_modres(res_exo, res_std, string, &n_exo);
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
      if(n_atom>0){
	// Write old residue
	res_i->n_atom=n_atom;
	// Write new residue
	nres++; res_i++;
      }
      if(nres >= MAXRES){
	printf("ERROR, too many residues (more than %d)\n", MAXRES);
	exit(8);
      }
      n_atom=0;
      if(res_i->atom==NULL)
	res_i->atom=malloc(MAXATOM*sizeof(struct atom));

      //res_i->aa=Get_aaname(string+17);
      char *s3=string+17;
      res_i->aa=Code_3_1(s3);
      if(res_i->aa=='X')res_i->aa=Het_res(s3, res_exo, res_std, n_exo);

      strcpy(res_i->label, res_label);
      // Store
      strcpy(res_old, res_label);
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

  for(int i=0; i<nres; i++){
    if(res[i].atom)free(res[i].atom); res[i].atom=NULL;
  }

  return(nres);
}

int Write_protein(struct protein *prot, struct residue *res, int nres,
		  char *filename, char chain_to_read, char exp_meth)
{
  struct atom *atom1, *atom2; 
  int i, j, nca=0, natoms=0;
  for(i=0; i<nres; i++)natoms+=res[i].n_atom;
  prot->natoms=natoms;

  //Set_prot_name(prot->code, filename, chain_to_read);
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
  prot->seqres=Read_seqres(prot, nres, filename, chain_to_read);


  return(0);
}

short *Read_seqres(struct protein *prot, int nres,
		   char *filename, char chain_to_read)
{
  FILE *file_in=fopen(filename, "r");
  if(file_in==NULL){
    printf("ERROR, file %s does not exist, exiting\n", filename);
    return(NULL);
  }
  printf("Reading seqres in PDB: %s\n",filename);
  short *seqres=malloc(nres*sizeof(short));
  char *seq3=NULL, string[300];
  int n_seqres=0, n=0, ndis=0;
  int dismax=nres, res_dis[dismax]; char aa_dis[dismax];

  while(fgets(string, sizeof(string), file_in)!=NULL){
    if(strncmp(string, "SEQRES", 6)==0){
      //SEQRES   1 A  334  MET GLY SER ASP LYS ILE HIS HIS HIS HIS HIS HIS GLU
      char chain=string[11];
      if(chain==chain_to_read){
	if(n_seqres==0){
	  sscanf(string+12, "%d", &n_seqres);
	  seq3=malloc(3*n_seqres*sizeof(char));
	}
	Get_seqres(string, seq3, &n);
      }
    }else if(strncmp(string, "REMARK 465", 10)==0){
      /* REMARK 465   M RES C SSSEQI
	 REMARK 465     MET A   -19 */
      if(strncmp(string+11, "IDENTI", 6)==0){continue;}
      if(strncmp(string+11, "EXPERI", 6)==0){continue;}
      if(strncmp(string+11, "THE FO", 6)==0){continue;}
      if(string[19]==chain_to_read && string[13]!='M'){
	if(ndis>=dismax){
	  printf("WARNING, too many disordered residues (>%d) Just read:\n",
		 dismax, string); continue;
	}
	sscanf(string+23, "%d", res_dis+ndis);
	printf("Dis. res. %s", string);
	aa_dis[ndis]=Code_3_1(string+15);
	ndis++;
      }
    }else if(strncmp(string, "ATOM", 4)==0){
       break;
    }
  }
  fclose(file_in);
  printf("%d SEQRES residues %d structured, %d disordered found in %s\n",
	 n_seqres, nres, ndis, filename);

  // The SEQRES record was not found, copy from ATOM
  int is_seqres=1, error=0;
  if(n_seqres==0){
    is_seqres=0;
    printf("WARNING, the SEQRES record is not present in %s\n", filename);
    printf("I will copy it from the ATOM record (%d res) omitting "
	   "%d disordered residues\n", nres, ndis);
    n_seqres=nres+ndis;
    n=n_seqres;
  }else if(n!=n_seqres){
    printf("ERROR reading SEQRES of chain %c, nres=%d expected %d\n",
	   chain_to_read, n, n_seqres); error=1;
  }

  // Get aa of SEQRES in 1 letter code
  char *seqr=malloc(n_seqres*sizeof(char)); int i;
  prot->seqr=seqr;
  for(i=0; i<n_seqres; i++){seqr[i]='X';}
  if(seq3){
    char *s1=seqr, *s3=seq3; int n_mod=0;
    for(i=0; i<n_seqres; i++){
      *s1=Code_3_1(s3);
      if(*s1=='X'){
	*s1=Het_res(s3, res_exo, res_std, n_exo); n_mod++;
      }
      printf("%c", *s1);
      s1++; s3+=3;
    }
    if(n_mod)printf(" %d modified residues", n_mod);
    printf("\n");
    int k=0, error=0;
    for(i=0; i<nres; i++){
      while(seqr[k]!=prot->aseq[i] && k<n_seqres){k++;}
      if(k>=n_seqres){
	printf("WARNING, i=%d a.a. %c not found\n", i, prot->aseq[i]);
	error++;
      }else{
	seqres[i]=k;
      }
    }
    if(error>10){
      printf("seqres: ");
      for(i=0; i<n_seqres; i++){printf("%c", seqr[i]);}
      printf("\natomres: ");
      for(i=0; i<nres; i++){printf("%c", prot->aseq[i]);}
      printf("\n");
      printf("ERROR, >= %d a.a. not found in seqres, exiting\n",error);
      exit(8);
    }
  }else{ // seqres was not found
    int k=0, j=-1, error=0; 
    for(i=0; i<nres; i++){
      int r=Get_number(prot->pdbres[i]);
      for(j=0; j<ndis; j++){
	if(res_dis[j]>=r){break;}
	seqr[k]=aa_dis[j]; k++;
      }
      if(k<n_seqres){
	seqr[k]=prot->aseq[i];
	seqres[i]=k; k++;
      }else{
	error++;
      }
    }
    if(error>10){
      printf("ERROR, too many seqres residues %d = %d+%d excess= %d\n",
	     n_seqres, nres, ndis, error); exit(8);
    }
  }
  return(seqres);
}

void Get_seqres(char *string, char *seq3, int *n_seqres)
{
  char *ptr=string+19; 
  /* If amino acid chain, aa=1 */
  char *seq=seq3+3*(*n_seqres);
  for(int k=0; k<13; k++){
    for(int j=0; j<3; j++){*seq=*ptr; seq++; ptr++;}
    (*n_seqres)++; ptr++; if(*ptr==' ')return;
  }
}

int Read_modres(char **res_exo, char **res_std, char *string, int *n_exo)
{
  for(int j=0; j<*n_exo; j++){
    if(strncmp(res_exo[j], string+12, 3)==0)return(0);
  }
  res_exo[*n_exo]=malloc(3*sizeof(char));
  res_std[*n_exo]=malloc(3*sizeof(char));
  res_exo[*n_exo][0]=string[12]; res_std[*n_exo][0]=string[24];
  res_exo[*n_exo][1]=string[13]; res_std[*n_exo][1]=string[25];
  res_exo[*n_exo][2]=string[14]; res_std[*n_exo][2]=string[26];

  (*n_exo)++; return(0);
} 


char Het_res(char *res_type_old, char **res_exo, char **res_std, int n_exo)
{
  for(int i=n_exo-1; i>=0; i--){
      if(strncmp(res_type_old,res_exo[i],3)==0){
	return(Code_3_1(res_std[i]));
      }
  }
  return('X');
}

char Code_3_1(char *res){

  if(strncmp(res,"ALA",3)==0){return('A');
  }else if(strncmp(res,"GLU",3)==0){return('E');
  }else if(strncmp(res,"GLN",3)==0){return('Q');
  }else if(strncmp(res,"ASP",3)==0){return('D');
  }else if(strncmp(res,"ASN",3)==0){return('N');
  }else if(strncmp(res,"LEU",3)==0){return('L');
  }else if(strncmp(res,"GLY",3)==0){return('G');
  }else if(strncmp(res,"LYS",3)==0){return('K');
  }else if(strncmp(res,"SER",3)==0){return('S');
  }else if(strncmp(res,"VAL",3)==0){return('V');
  }else if(strncmp(res,"ARG",3)==0){return('R');
  }else if(strncmp(res,"THR",3)==0){return('T');
  }else if(strncmp(res,"PRO",3)==0){return('P');
  }else if(strncmp(res,"ILE",3)==0){return('I');
  }else if(strncmp(res,"MET",3)==0){return('M');
  }else if(strncmp(res,"PHE",3)==0){return('F');
  }else if(strncmp(res,"TYR",3)==0){return('Y');
  }else if(strncmp(res,"CYS",3)==0){return('C');
  }else if(strncmp(res,"TRP",3)==0){return('W');
  }else if(strncmp(res,"HIS",3)==0){return('H');
  }else if(strncmp(res,"HIE",3)==0){return('H');
  }else if(strncmp(res,"HID",3)==0){return('H');
  }else if(strncmp(res,"HIP",3)==0){return('H');
  }else if(strncmp(res,"ASX",3)==0){return('N');
  }else if(strncmp(res,"GLX",3)==0){return('Q');
  }else{
    char MODRES=Get_aaname(res);
    if(MODRES=='X'){
      MODRES=Het_res(res, res_exo, res_std, n_exo);
    }
    if(MODRES=='X'){
      printf("WARNING, a.a. %c%c%c not known\n", *res,*(res+1),*(res+2));
      //exit(8);
    }
    return(MODRES);
  }
}

int Get_number(char *pdbres){
  // Find number associated to pdbres
  int r, ic=0;
  char res[10], *cr=res; strcpy(res, pdbres);
  while(*cr!='\0' && ic<9){
    if(isdigit(*cr)==0){*cr=' '; break;} cr++; ic++;
  }
  *cr=' ';
  sscanf(res, "%d", &r);
  return(r);
}

int Select_domain(struct protein *prot, 
		  int *ini_frag, int *end_frag, int nfrag)
{
  /* The numbering of cath domains seem to be based on SEQRES, not on
     the RES record in PDB nor on structured residues.
     Ex: 4zi8 A 227-334 SEQRES: 1-334 RES: -19-314 Str: 2-314  */
  int Lsel=0;
  printf("Selecting %d fragments:", nfrag);
  for(int i=0; i<nfrag; i++){
    printf(" %d-%d", ini_frag[i], end_frag[i]);
    Lsel+=end_frag[i]-ini_frag[i]+1;
  }
  printf(" Total: %d\n", Lsel);

  int L=prot->len, select[L], i, sel=0;
  for(i=0; i<L; i++)select[i]=0;

  for(i=0; i<L; i++){
    int r=prot->seqres[i]; //i+1;
    if(0){r=Get_number(prot->pdbres[i]);}
    for(int k=0; k<nfrag; k++){
      if(r>=ini_frag[k] && r<=end_frag[k]){select[i]=1; sel++; break;}
    }
  }
  //printf("Selecting %d residues out of %d in %d fragments\n", sel,L,nfrag);
    int tol=Lsel/3;
    if(sel<Lsel){
    printf("WARNING, %d residues found in domain, expected %d\n", sel, Lsel);

    if(sel<Lsel-tol && sel<L-tol){
      printf("seqres: ");
      for(i=0; i<L; i++){printf("%d ", prot->seqres[i]);}
      printf("\nEXITING\n"); exit(8);
    }
  }

  sel=0;
  float *xca_rot=prot->xca_rot;
  for(i=0; i<L; i++){
    if(select[i]==0){continue;}
    if(i==sel){sel++; continue;}
    // copy i to sel
    float *xca_old=prot->xca_rot+3*i;
    *(xca_rot)=*(xca_old); xca_rot++; xca_old++;
    *(xca_rot)=*(xca_old); xca_rot++; xca_old++;
    *(xca_rot)=*(xca_old); xca_rot++; xca_old++;
    prot->aseq[sel]=prot->aseq[i];
    prot->n_atom[sel]=prot->n_atom[i];
    strcpy(prot->pdbres[sel], prot->pdbres[i]);
    prot->seqres[sel]=prot->seqres[i];
    if(prot->seq1)prot->seq1[sel]=prot->seq1[i];
    if(prot->seq_bs)prot->seq_bs[sel]=prot->seq_bs[i];
    prot->res_atom[sel]=prot->res_atom[i];
    int s3=3*sel, i3=3*i;
    if(prot->vec){
      prot->vec[s3]=prot->vec[i3]; s3++; i3++;
      prot->vec[s3]=prot->vec[i3]; s3++; i3++;
      prot->vec[s3]=prot->vec[i3];
    }
    sel++;
  }
  printf("Selecting %d residues out of %d in %d fragments\n", sel,L,nfrag);
  prot->nca=sel;
  prot->len=sel;
  return(sel);
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
  if(chain != tmp[4]){
    if(i<CHARPDB){if(chain!=' '){tmp[i]=chain;}else{tmp[i]='_';} i++;}
  }
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
  if(N_sec_ele){
    for(i_res=0; i_res<N_res; i_res++)printf("%c",SS_code[ss3[i_res]]);
    printf("\n");
  }
  //exit(8);
  return(N_sec_ele);
}

void Modres(){
  int i;
  for(i=0; i<20; i++)modres[i]=malloc(200*sizeof(char));
  strcpy(modres[0], "CSDMDOORNDBZNALLALHACAYANCB"); // Ala
  strcpy(modres[1], // Cys 
	 "CYGCY3CSUCSPCYMCASCSBCSRCMECMHCSSCSXCSWCSOALSSMCCEAOCSCCSSCYYCM");
  strcpy(modres[2], "BFDSNNASABHDSUI"); // Asp
  strcpy(modres[3], "PCAAR4GLQGMA");    // Glu 
  strcpy(modres[4], "DPN"); // Phe
  strcpy(modres[5], "GL3CHGGLZACYFGL"); //  Gly 
  strcpy(modres[6], "MHSNEPHIC"); //  His
  strcpy(modres[7], "IIL"); // Ile
  strcpy(modres[8], "LLPKCXLYZALYLCXMCLMLYMLZM3L"); // Lys 
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
  printf("WARNING, a.a. %c not known\n", aseq);
  return(20);
}


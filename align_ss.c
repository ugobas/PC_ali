#include "protein.h"
#include "align_ss.h"
#include "Contact_divergence_aux.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <math.h>


int *warn_notali=NULL;
int verbose=1;
extern char SS_code[];

static int SS_match(int **Ali, short **ss_prot, int m, int i,
		    int ini_res, int end_res, int move_res);
static int Copy_ali_nogaps(int **column, int **Ali_new,
			   int **Ali_store, int *L, int m, int N_ali);
static int gap_score(int *Ali_i, int *col_gaps,
		     int ini_gap, int end_gap, int ini_new);
static void Count_gaps(int *col_gaps, int **msa, int m, int N_ali);


int Set_sec_str(struct protein *prots, int N_pdb)
{
  // Determine number of secondary structures
  float ss_ave=0;
  for(int i=0; i<N_pdb; i++){
    if(prots[i].ss3==NULL){
      printf("ERROR, secondary structure of %s (i=%d) not present\n",
	     prots[i].domname, i); return(-1);
    }
    short *ss=prots[i].ss3; int nss=0; short type=0;
    for(int j=0; j<prots[i].len; j++){
      if(ss[j]!=0 && ss[j]!=type){nss++; type=ss[j];}
    }
    ss_ave+=nss;
  }
  printf("Mean number of secondary structure elements: %.1f\n",ss_ave/N_pdb);
  return(0);
}


void Write_ali_prot(int **msa, struct protein **prots,
		    int N_pdb, int N_ali, char *name_in, char *what)
{
  char f_ali_ss[100];
  sprintf(f_ali_ss,"%s_%s.msa",name_in, what);
  printf("Writing msa in %s\n", f_ali_ss);
  FILE *f_ali=fopen(f_ali_ss,"w");
  for(int i=0; i<N_pdb; i++){
    struct protein *prot=prots[i];
    fprintf(f_ali,">%s\n",prot->domname);
    for(int j=0; j<N_ali; j++){
      if(msa[i][j]<0){
	fprintf(f_ali,"-");
      }else{
	fprintf(f_ali,"%c",prot->aseq[msa[i][j]]);
      }
    }
    fprintf(f_ali,"\n");
  }
  fclose(f_ali);
}

void Write_ss_ali(int **ali, struct protein **prot_p, int N_pdb, int N_ali,
		  char *name, char *what)
{
  char f_ali_ss[100]; sprintf(f_ali_ss,"%s.%s_ss.msa",name,what);
  printf("Writing secondary structure MSA in %s\n", f_ali_ss);
  FILE *f_ali=fopen(f_ali_ss,"w");
  for(int i=0; i<N_pdb; i++){
    fprintf(f_ali,">%s\n",prot_p[i]->domname);
    short *ss3=prot_p[i]->ss3;
    for(int j=0; j<N_ali; j++){
      if(ali[i][j]<0){
	fprintf(f_ali,"-");
      }else{
	fprintf(f_ali,"%c",SS_code[ss3[ali[i][j]]]);
      }
    }
    fprintf(f_ali,"\n");
  }
  fclose(f_ali);

}

int Align_ss(int **Ali_new,
	     int **Ali_store, int N_ali, struct protein **prots, int m,
	     int SHIFT_MAX, int score_type)
{

  /* For every column in MSA Ali_store, test whether there is a gap in a seq and a residue in the core of a secondary str. on another one. In this case, determine the initial and final column of the secondary structure. Compute the score of moving the gaps present in these columns that do not coincide with the tails of the secondary stucture either to the initial or to the final column. Perform the change with minimum score. Move to the final column of the secondary structure. */

  if(verbose)
    printf("Modifying MSA by eliminating gaps in secondary structures\n");

  // Store Alignments
  int L[m], Test_gap[m], *column[m], i; short *ss3_prot[m];
  for(i=0; i<m; i++){
    column[i]=malloc(N_ali*sizeof(int));
    ss3_prot[i]=(prots[i])->ss3;
    L[i]=(prots[i])->len;
    Test_gap[i]=1;
  }
  int N_ali_nogap=Copy_ali_nogaps(column, Ali_new, Ali_store, L, m, N_ali);
  int col_gaps[N_ali_nogap];
  Count_gaps(col_gaps, Ali_new, m, N_ali_nogap);

  int nmove=0;
  for(int icol=0; icol<N_ali_nogap; icol++){
    // Is there a gap?
    if(col_gaps[icol]==0)continue;
    int gap=0;
    for(i=0; i<m; i++){
      if(Ali_new[i][icol]<0){if(Test_gap[i])gap=1;}
      else{Test_gap[i]=1;}
    }
    if(gap==0)continue;
    // Yes, there is a gap at icol. Is it inside a secondary structure?
    int ini_ss=-1, end_ss=-1;
    int ini_ssi[m], end_ssi[m];
    for(i=0; i<m; i++){
      ini_ssi[i]=-1; end_ssi[i]=-1;
      int iss=Ali_new[i][icol], j=icol;
      if(iss<0){
	while(j<N_ali_nogap){
	  if(Ali_new[i][j]>=0){iss=Ali_new[i][j]; break;} j++;
	}
	if(iss<0)continue;
      }
      short *ss=ss3_prot[i]+iss; 
      if((*ss!=1)&&(*ss!=2))continue; // It is a loop
      if((j>icol)&&(*ss!=*(ss-1)))continue; // leave if gap and no ss.
      int *col=column[i], iss2=iss; short *ss2=ss;
      ini_ssi[i]=col[iss];
      while((iss2>=0)&&(*ss2==*ss)){
	if((col[iss2]<ini_ssi[i])&&(col[iss2]>=0))ini_ssi[i]=col[iss2];
	iss2--; ss2--;
      }
      iss2=iss; ss2=ss;
      end_ssi[i]=col[iss];
      while((iss2<L[i])&&(*ss2==*ss)){
	if(col[iss2]>end_ssi[i])end_ssi[i]=col[iss2];
	iss2++; ss2++;
      }
      if((ini_ss<0)||(ini_ssi[i]<ini_ss))ini_ss=ini_ssi[i];
      if((end_ss<0)||(end_ssi[i]>end_ss))end_ss=end_ssi[i];
    }
    if((ini_ss<0)||(end_ss<0)){ //||(ini_ss==icol)||(end_ss==icol)){
      // Accept gap at extreme of secondary structure element
      for(i=0; i<m; i++){if(Ali_new[i][icol]<0)Test_gap[i]=0;}
      continue;
    }

    int last_col=end_ss; int idir_old=0;
    for(i=0; i<m; i++){
      if(Test_gap[i]==0)continue;
      // Gap inside secondary structure, move toward ini_ss or end_ss
      // Identify all gaps and initial and final residues
      if(0){
	// Find boundaries of secondary structure excluding target seq
	ini_ss=-1; end_ss=-1;
	for(int j=0; j<m; j++){
	  if(j==i)continue;
	  if((ini_ss<0)||(ini_ssi[j]<ini_ss))ini_ss=ini_ssi[j];
	  if((end_ss<0)||(end_ssi[j]>end_ss))end_ss=end_ssi[j];
	}
	if(ini_ss<0)continue;
      }

      int *col=column[i], *ali_i=Ali_new[i], j, jj;
      int ngap=0, open=0, inigap[N_ali], endgap[N_ali];
      int ini_res=-1, end_res=-1, next_ss=-1;
      short *ss_ini=NULL, *ss=ss3_prot[i];
      for(j=icol; j<=end_ss; j++){
	if(ali_i[j]<0){
	  if(open==0){inigap[ngap]=j; open=1;}
	  endgap[ngap]=j;
	}else{
	  if(ini_res<0)ini_res=j;
	  end_res=j;
	  if(open){open=0; if(1 || inigap[ngap]!=ini_ss)ngap++;}
	  short *ssj=ss+ali_i[j];
	  if(ss_ini==NULL){ss_ini=ssj;}
	  else if((*ss_ini!='H')&&(*ss_ini!='E')){ss_ini=ssj;}
	  else if((*ssj!=*ss_ini)&&((*ssj=='H')||(*ssj=='E')))
	  {next_ss=j;break;}
	}
      }
      if(open){open=0; if(1 || endgap[ngap]!=end_ss)ngap++;}
      if(0){
	printf("col %d seq %d %d gaps:", icol, i, ngap);
	for(j=0; j<ngap; j++){printf(" %d-%d", inigap[j], endgap[j]);}
	printf("\n");
      }

      if((next_ss>=0)&&(next_ss<last_col))last_col=next_ss;
      // Backward
      if(ali_i[icol]<0){
	j=icol-1; while(ali_i[j]<0){inigap[0]=j; j--;}
      }
      for(j=icol-1; j>=ini_ss; j--){
	if(ali_i[j]>=0){
	  short *ssj=ss+ali_i[j];
	  if(ss_ini==NULL){ss_ini=ssj;}
	  else if((*ss_ini!='H')&&(*ss_ini!='E')){ss_ini=ssj;}
	  else if((*ssj!=*ss_ini)&&((*ssj=='H')||(*ssj=='E'))){break;}
	  ini_res=j;
	}
      }
      if((ngap==0)||(ini_res<0)){Test_gap[i]=0; continue;}
 
      // Sequence i has ngap gaps inside secondary structure
      // Determine shifts to the left and to the right
      int lgap[ngap];
      for(j=0; j<ngap; j++)lgap[j]=endgap[j]-inigap[j]+1;

      int store_l=ini_res, store_r=end_res;
      char dir[10]; int idir=0;
 
      // recursively move the gap with lowest shift
      int todo[ngap]; for(j=0; j<ngap; j++)todo[j]=1;
      for(jj=0; jj<ngap; jj++){

	int inimin=1000, igap1=-1, endmax=-1, igap2=-1;
	for(j=0; j<ngap; j++){
	  if(todo[j]==0)continue;
	  if(igap1<0 || inigap[j]<inimin){inimin=inigap[j]; igap1=j;}
	  if(igap2<0 || endgap[j]>endmax){endmax=endgap[j]; igap2=j;}
	}
	int igap=-1, res1=-1, res2=-1;
	int shift_min=1000, score_max=-10, score=0, store_res=0, gap1=0;
	if(igap1>=0){
	  /* gap is moved to the left between inires_mv and inires_mv+lgap.
	       Residues between inires_mv and inigap-1 are moved to the right
	       between inires_mv+lgap and endgap */
	  int endres_mv=inigap[igap1]-1;
	  for(int inires_mv=store_l; inires_mv>=0; inires_mv--){
	    if(ali_i[inires_mv]<0)break;
	    int inires_new=inires_mv+lgap[igap1]; // Where residues are moved
	    if(score_type){
	      score=SS_match(Ali_new,ss3_prot,m,i,
			     inires_mv, endres_mv, inires_new);
	    }else{
	      score=gap_score(Ali_new[i],col_gaps,
			      inigap[igap1],endgap[igap1],inires_mv);
	    }
	    if((igap)&&(score<score_max))break;
	    if((igap<0)||(score>score_max)||
	      (score==score_max && idir_old<0)){
	      igap=igap1; strcpy(dir,"left"); idir=-1;
	      score_max=score;
	      store_res=inires_new;
	      res1=ali_i[inires_mv];
	      res2=ali_i[endres_mv];
	      gap1=inires_mv; // Where gap is moved
	      shift_min=inigap[igap]-gap1;
	    }
	  }
	}
	if(igap2>=0){
	    /* gap is moved to the right between inires_mv and inires_mv+lgap.
	       Residues between inires_mv and inigap-1 are moved to the right
	       between inires_mv+lgap and endgap */
	  int score_right=-1;
	  int inires_new=inigap[igap2], inires_mv=endgap[igap2]+1;
	  for(int endres_mv=store_r; endres_mv<N_ali_nogap; endres_mv++){
	    if(ali_i[endres_mv]<0)break;
	    if(score_type){
	      score=SS_match(Ali_new,ss3_prot,m,i,
			     inires_mv, endres_mv, inires_new);
	    }else{
	      score=
		gap_score(Ali_new[i],col_gaps,
			  inigap[igap2],endgap[igap2],endres_mv-lgap[igap2]+1);
	    }
	    //if((igap<0)||(shift<shift_min)){
	    if((score_right>=0)&&(score<score_right))break;
	    score_right=score;
	    if((igap<0)||(score > score_max)||
	      (score==score_max && idir_old>0)){
	      igap=igap2; strcpy(dir,"right"); idir=1;
	      score_max=score;    
	      store_res=inires_new;
	      res1=ali_i[inires_mv];
	      res2=ali_i[endres_mv];
	      gap1=endres_mv-lgap[igap]+1;
	      shift_min=gap1-inigap[igap];
	    }
	  }
	}
	if(igap<0)continue;
	if(shift_min<=0){todo[igap]=0; continue;}
	if(shift_min >SHIFT_MAX){
	  if(verbose)
	    printf("WARNING, icol %d seq %d gap %d of %d shift=%d > %d, "
		   "keeping ali\n",
		   icol,i,jj+1, ngap, shift_min, SHIFT_MAX);
	  //printf("gap: %d-%d store_l= %d store_r=%d\n",
	  //	 inigap[igap], endgap[igap], store_l, store_r);
	  break;
	}
	todo[igap]=0;
	idir_old=idir;

	if(res1<0 || res2<0 || (res2<res1)){
	  printf("WARNING gap %d to %s, res1=%d res2=%d exiting\n",
		 jj+1,dir,res1,res2);
	  for(j=0; j<ngap; j++){
	    printf("gap %d: %d-%d\n",j,inigap[j],endgap[j]);
	  }
	  continue;
	}

	// The move is accepted
	if(strcmp(dir,"left")==0){
	  // Shift to the left
	  store_l=store_res;
	}else{
	  store_r=gap1-1;
	}
	// Move residues and gap
	for(int k=res1; k<=res2; k++){
	  col[k]=store_res; ali_i[store_res]=k; store_res++;
	}
	for(j=gap1; j<(gap1+lgap[igap]); j++){ali_i[j]=-1; col_gaps[j]++;}
	for(j=inigap[igap]; j<=endgap[igap]; j++){col_gaps[j]--;}
	if(verbose){
	  printf("icol %d seq %d "
		 "moving gap %d l=%d from col %d-%d to %d-%d shift= %d\n",
		 icol,i,igap+1,lgap[igap],inigap[igap],endgap[igap],
		 gap1,gap1+lgap[igap]-1,idir*shift_min);
	}	
	nmove++;

	// Test and print

	/*if(shift_r[igap]<shift_l[igap]){
	  // Test (right)
	  k=store_r-lgap[igap]; int ak=ali_i[k];
	  if((ak>=0)&&(ak<(N_ali_nogap-1))&&(col[ak]>=col[ak+1])){
	    printf("ERROR moving to %s (end), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak+1, col[ak+1]); exit(8);
	  }
	  k=inigap[igap]; ak=ali_i[k];
	  if((ak>0)&&(col[ak]<=col[ak-1])){
	    printf("ERROR moving to %s (ini), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak-1, col[ak-1]); exit(8);
	  }
	}else{
	  // Test left
	  k=store_l+lgap[igap]; int ak=ali_i[k];
	  if((ak>0)&&(col[ak]<=col[ak-1])){
	    printf("ERROR moving to %s (ini), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak-1, col[ak-1]); exit(8);
	  }
	  k=endgap[igap]; ak=ali[k];
	  if((ak>=0)&&(ak<(N_ali_nogap-1))&&(col[ak]>=col[ak+1])){
	    printf("ERROR moving to %s (end), col[%d]=%d col[%d]=%d\n",
		   dir, ak, col[ak], ak+1, col[ak+1]); exit(8);
	  }
	  }*/
	int ires=ali_i[ini_res];
	if((ires>0) && (col[ires]<=col[ires-1])){
	  printf("ERROR moving gaps to %s  col[%d]=%d col[%d]=%d\n",
		 dir, ires, col[ires], ires-1, col[ires-1]);
	  printf("sequence= %d gap %d of %d l=%d shift= %d\n",
		 i, jj, ngap, lgap[igap], shift_min);
	  printf("col %d columns of sec.str: %d-%d sec.str. residues: %d-%d\n",
		 icol, ini_ss, end_ss, ires, ali_i[end_res]);
	  printf("Columns of residues: ");
	  for(j=ires; j<=ali_i[end_res]; j++)printf(" %d", col[j]);
	  printf("\n");
	  if(ini_res)printf("Previous column of residues: %d\n",col[ires-1]);
	  exit(8);
	}
      } // end gaps
    } // end sequences
    if(last_col>icol){icol=last_col; printf("icol= %d)\n", icol);}
    // All columns until last_col have been changed
  }

 for(i=0; i<m; i++){
   int *col=column[i], j;
    for(j=0; j<L[i]; j++){
      if((col[j]>=N_ali_nogap)){
	printf("ERROR, sequence %d res %d column %d not in [0,%d]\n",
	       i, j, col[j], N_ali_nogap-1); exit(8);
      }else if(j && (col[j]<=col[j-1])&&(col[j-1]>=0)&&(col[j]>=0)){
	printf("ERROR in seq %d, column[%d]=%d<=column[%d]=%d\n",
	       i, j, col[j], j-1, col[j-1]); exit(8);
      }

    }
  }
 int N_new_nogap=  
   Copy_ali_nogaps(column, Ali_new, Ali_new, L, m, N_ali_nogap);

 if(verbose){
   printf("%d gaps moved, length of alignment changed from %d to %d\n",
	  nmove, N_ali_nogap, N_new_nogap);
 }

 // Free memory
 for(i=0; i<m; i++){free(column[i]);}
 //exit(8);

 return(N_new_nogap);
}

int SS_match(int **Ali, short **ss_prot, int m, int i,
	     int ini_res, int end_res, int new_res)
{
  // Fragment of prot i moved from columns ini_res-end_res minus
  // to new_res-new_res+l_res. Compute new match minus previous one
  int ss_match=0; int new=new_res;
  int *Ali_i=Ali[i]; short *ss_i=ss_prot[i];
  for(int j=ini_res; j<=end_res; j++){
    int res=Ali_i[j]; if(res<0)continue;
    //if(ss_i[res]=='c')continue;
    for(int k=0; k<m; k++){
      if(k==i)continue;
      if((Ali[k][j]>=0)&&(ss_prot[k][Ali[k][j]]==ss_i[res]))ss_match--;
      if((Ali[k][new]>=0)&&(ss_prot[k][Ali[k][new]]==ss_i[res]))ss_match++;
    }
    new++;
  }
  return(ss_match);
}

int gap_score(int *Ali_i, int *col_gaps,
	      int ini_gap, int end_gap, int ini_new)
{
  // move gap of seq.i from ini_gap to ini_new.
  // Count max. value of col_gaps
  int max_old=-1, max_new=-1;
  int n=ini_new;
  for(int j=ini_gap; j<=end_gap; j++){
    int col_j=col_gaps[j], col_n=col_gaps[n];
    if(col_j>max_old)max_old=col_j;
    if(col_n>max_old)max_old=col_n;
    if(Ali_i[j]<0){col_j--; col_n++;}
    if(Ali_i[n]<0){col_n--; col_j++;}
    if(col_j>max_new)max_new=col_j;
    if(col_n>max_new)max_new=col_n;
    n++;
  }
  return(max_new-max_old);
}

int Test_notali(int **Ali, struct protein *prot, int n, int N_ali)
{
   // Test if some structural position is not aligned
  int n_not=0, s_not=0;
  int col[N_ali], i, j;
  for(i=0; i<n; i++){
    int L=prot[i].len;
    for(j=0; j<N_ali; j++)col[j]=-1;
    for(j=0; j<N_ali; j++)if(Ali[i][j]>=0)col[Ali[i][j]]=j;
    int not_i=0; int notali[L];
    for(j=0; j<L; j++)if(col[j]<0){notali[not_i]=j; not_i++;}
    if(not_i){
      printf("WARNING, protein %d L=%d %d columns not aligned: ",i,L,not_i);
      for(j=0; j<not_i; j++)printf(" %d",notali[j]);
      printf("\n"); n_not+=not_i; s_not++;
    }
  }
  if(n_not){
    printf("WARNING, %d residues not aligned in %d proteins\n",n_not, s_not);
    if(0){printf("Exiting\n"); exit(8);}
  }
  return(s_not);
}

int Copy_ali_nogaps(int **column, int **Ali_new,
		    int **Ali_store, int *L, int m, int N_ali)
 {
   int N_ali_nogap=N_ali,  all_gap[N_ali], i, j;
   for(j=0; j<N_ali; j++){
     all_gap[j]=1;
     for(i=0; i<m; i++){
      if(Ali_store[i][j]>=0){all_gap[j]=0; break;}
    }
    if(all_gap[j])N_ali_nogap--;
  }
  for(i=0; i<m; i++){
    int *col=column[i]; for(j=0; j<N_ali; j++){col[j]=-1;}
    int nali=0;
    int *al=Ali_store[i];
    int *ali=Ali_new[i];
    int len=L[i];

    for(j=0; j<N_ali; j++){
      if(all_gap[j])continue;
      ali[nali]=al[j];
      if(ali[nali]>len)ali[nali]=-1;
      if(ali[nali]>=0)col[ali[nali]]=nali;
      nali++;
    }
    for(j=nali; j<N_ali; j++){ali[j]=-1;}
  }
  return(N_ali_nogap);
}

void Count_gaps(int *col_gaps, int **msa, int m, int N_ali){
  int j; for(j=0; j<N_ali; j++)col_gaps[j]=0;
  for(int i=0; i<m; i++){
    int *ali=msa[i];
    for(j=0; j<N_ali; j++){if(*ali<0)col_gaps[j]++; ali++;}
  }
}

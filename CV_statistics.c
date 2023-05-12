#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "CV_statistics.h"

int Ndist=3;
float *alpha_fit, *thr_fit;

static void Ave_se(double *x1, double *x2, int n);
static void Print_line(char *out, float *alpha_val, int Na,
		       double *Acc1, double *Acc2, int *num_signi,
		       int num, double *xd1, double *xd2,
		       double *t1, double *t2,
		       double *TV1, double *Sign1, double *Out1);
float Fit_alpha(float *thr,
		float **Div, float **CV, float **t, int **nout,
		float **Seq_Id, float SI_thr, int N,
		FILE *file_out);

void CV_statistics(char *name_in, char *head, char *OUTG, char *ALI,
		   char *name_div, float **div, float **CV, float **t,
		   int **nout, int **nsign, int **nTIV,
		   float **Seq_Id, float SI_thr,
		   float **fun_sim, float fun_low, float fun_high,
		   int N, int dist, char *funtype)
{
  // Parameters
  int Nalpha=6, Na1=Nalpha-1, Na=Na1; // Number of alpha values
  float alpha_val[Nalpha], Acc_thr[Nalpha];
  for(int a=0; a<Nalpha; a++)Acc_thr[a]=0.15;
  alpha_val[0]=0.5; alpha_val[1]=0.7; alpha_val[2]=1.0;
  //alpha_val[0]=0.0; alpha_val[1]=0.5; alpha_val[2]=1.0;
  if(Nalpha>4)alpha_val[3]=1.5;
  if(Nalpha>5)alpha_val[4]=2.0;

  // Make and print fit
  char name_out[100]; FILE *file_out=NULL;
  int FIT=0;
  if(FIT){
    if(fun_sim==NULL){
      sprintf(name_out, "%s_%s_Outg%s_fit_alpha.dat", name_in, ALI, OUTG);
      if(dist==0){
	file_out=fopen(name_out, "w");
	printf("Writing %s\n", name_out);
	alpha_fit=malloc(Ndist*sizeof(float));
	thr_fit=malloc(Ndist*sizeof(float));
      }else{
	file_out=fopen(name_out, "a");
      }
      fprintf(file_out, "# Divergenge: %s\n", name_div);
      if(dist<Ndist){
	alpha_fit[dist]=
	  Fit_alpha(thr_fit+dist, div, CV, t, nout, Seq_Id, SI_thr, N, file_out);
      }else{
	printf("WARNING, not allowed distance %s with index=%d < %d\n",
	       name_div, dist, Ndist);
	printf("Increase Ndist in CV_statistics.c\n");
	fprintf(file_out, "WARNING, not allowed distance %s with index=%d < %d\n",
		name_div, dist, Ndist);
	fprintf(file_out, "Increase Ndist in CV_statistics.c\n");
	float dumm;
	Fit_alpha(&dumm, div, CV, t, nout, Seq_Id, SI_thr, N, file_out);
      }
      fclose(file_out);
    }
    if(dist<Ndist){
      alpha_val[Na1]=alpha_fit[dist];
      Acc_thr[Na1]=thr_fit[dist];
      Na=Nalpha;
    }
  }// End fit
      
  //float t_thr=3;
  int Numd=300;     // Maximum number of bins of distances
  float maxd=3;     //Maximum distance
  float d_step=maxd/Numd;
  int min_num=90;   // Minimum number of pairs in a bin of distances

  // Prepare counters
  int num_d[Numd], k, ia;
  float Acc_thr2[Na];
  double Acc_d1[Na][Numd], Acc_d2[Na][Numd], Significant_d[Na][Numd];
  double d1[Numd], d2[Numd], t_d1[Numd], t_d2[Numd];
  double Sign[Numd], Out[Numd], TV_d1[Numd];
  for(ia=0; ia<Na; ia++)Acc_thr2[ia]=Acc_thr[ia]*Acc_thr[ia];
  for(k=0; k<Numd; k++){
    num_d[k]=0; d1[k]=0; d2[k]=0;
    t_d1[k]=0; t_d2[k]=0; TV_d1[k]=0; Sign[k]=0; Out[k]=0;
    for(ia=0; ia<Na; ia++){
      Acc_d1[ia][k]=0; Acc_d2[ia][k]=0; Significant_d[ia][k]=0;
    }
  }
  // Loop over all pairs
  int discard=0;
  for(int i=0; i<N; i++){
    for(int j=0; j<i; j++){
      if(fun_sim && (fun_sim[i][j]>=0) &&
	 ((fun_sim[i][j]<fun_low)||(fun_sim[i][j]>fun_high)))
	continue;
      if(nout[i][j]==0){discard++; continue;}
      float d=div[i][j];
      float CVij=fabs(CV[i][j]);
      float tt=fabs(t[i][j]);
      int k=d/d_step; if(k>=Numd)k=Numd-1;
      num_d[k]++;
      for(ia=0; ia<Na; ia++){
	float Acc=CVij/pow(d,alpha_val[ia]);
	float Acc2=Acc*Acc;
	Acc_d1[ia][k]+=Acc;
	Acc_d2[ia][k]+=Acc2;
	if(Acc2*(1-1./(tt*tt))>Acc_thr2[ia]){
	  Significant_d[ia][k]++;
	}
      } // end alpha
      d1[k]+=d;
      d2[k]+=d*d;
      t_d1[k]+=tt;
      t_d2[k]+=tt*tt;
      TV_d1[k]+=nTIV[i][j];
      Sign[k]+=nsign[i][j];
      Out[k]+=nout[i][j];
    } //end j
  } // end i

  // Compute global averages
  double Acc1[Na], Acc2[Na]; int num_signi[Na], num=0;
  double xd1=0, xd2=0, t1=0, t2=0, TV1=0, Out1=0, Sign1=0;
  for(ia=0; ia<Na; ia++){
    Acc1[ia]=0; Acc2[ia]=0; num_signi[ia]=0;
  }
  for(k=0; k<Numd; k++){
    for(ia=0; ia<Na; ia++){
      Acc1[ia]+=Acc_d1[ia][k];
      Acc2[ia]+=Acc_d2[ia][k];
      num_signi[ia]+=Significant_d[ia][k];
    }
    num+=num_d[k];
    xd1+=d1[k];  xd2+=d2[k]; 
    t1+=t_d1[k]; t2+=t_d2[k];
    TV1+=TV_d1[k]; Sign1+=Sign[k];
    Out1+=Out[k];
  }
  char header[400], out[400], tmp[80];
  sprintf(header, "#1=<d> s.e."); k=3;
  for(ia=0; ia<Na1; ia++){
    sprintf(tmp, "\t%d=<acc>_%.2f s.e.\t%d=<signi> s.e.",
	    k, alpha_val[ia], k+2); k+=4;
    strcat(header, tmp);
  }
  if(Na>Na1){
    sprintf(tmp, "\t%d=<acc>_alpha s.e.\t%d=<signi> s.e.",
	    k, k+2); k+=4;
    strcat(header, tmp);
  }
  sprintf(tmp, "\t%d=t s.e.\t%d=Sign s.e.\t%d=TIV s.e.\t%d=Outg\t%d=num",
	  k, k+2, k+4, k+6, k+7);
  strcat(header, tmp);

  // Print global averages and print
  if(fun_sim==NULL){
    sprintf(name_out, "%s_%s_Outg%s_All_dist_%s.dat",
	    name_in, ALI, OUTG, funtype);
  }else{
    sprintf(name_out, "%s_%s_Outg%s_dist_%s_Functions.dat",
	    name_in, ALI, OUTG, name_div);
  }
  if(((fun_sim==NULL)&&(dist==0))||
     (fun_sim && (strcmp(funtype,"SameFun")==0))){
    file_out=fopen(name_out, "w");
    printf("Writing %s\n", name_out);
    fprintf(file_out, "%s", head);
    fprintf(file_out, "# %d pairs discarded\n", discard);
    fprintf(file_out, "%s\tdiv\n", header);
  }else{
    file_out=fopen(name_out, "a");
  }
  if(fun_sim)
    fprintf(file_out, "# Function similarity between %.3f and %.3f\n",
	    fun_low, fun_high);
  Print_line(out, alpha_val, Na, Acc1, Acc2, num_signi, 
	     num, &xd1, &xd2, &t1, &t2, &TV1, &Sign1, &Out1);
  fprintf(file_out, "%s\t%s\n", out, name_div);
  if(FIT)fprintf(file_out, "# alpha= %.2f Acc_thr= %.2f\n",
		 alpha_val[Na1], Acc_thr[Na1]);
  fclose(file_out);


  // Compute averages versus distance and print.
  if(fun_sim)return; // Do not distinguish functional classes
  sprintf(name_out, "%s_dist_%s_%s.dat", name_in, name_div, funtype);
  file_out=fopen(name_out, "w");
  printf("Writing %s\n", name_out);
  fprintf(file_out, "# divergence= %s %d pairs\n", name_div, num);
  fprintf(file_out, "%s\n", header);
  int ini=1;
  for(k=0; k<Numd; k++){
    if(ini){
      // Reset
      for(ia=0; ia<Na; ia++){
	Acc1[ia]=0; Acc2[ia]=0; num_signi[ia]=0;
      }
      num=0; xd1=0; xd2=0; t1=0; t2=0; TV1=0; Out1=0; Sign1=0;
    }
    // Sum
    num+=num_d[k];
    for(ia=0; ia<Na; ia++){
      Acc1[ia]+=Acc_d1[ia][k];
      Acc2[ia]+=Acc_d2[ia][k];
      num_signi[ia]+=Significant_d[ia][k];
    }
    xd1+=d1[k]; xd2+=d2[k];
    t1+=t_d1[k]; t2+=t_d2[k];
    TV1+=TV_d1[k]; Sign1+=Sign[k];
    Out1+=Out[k];
    if(num<min_num){ini=0; continue;}
    ini=1;
    Print_line(out, alpha_val, Na, Acc1, Acc2, num_signi, 
	       num, &xd1, &xd2, &t1, &t2, &TV1, &Sign1, &Out1);
    fprintf(file_out, "%s\n", out);
  }
  fclose(file_out);

}

void Print_line(char *out, float *alpha_val, int Na,
		double *Acc1, double *Acc2, int *num_signi,
		int num, double *xd1, double *xd2,
		double *t1, double *t2,
		double *TV1, double *Sign1,
		double *Out1)
{
  float p; char tmp[80];
  Ave_se(xd1, xd2, num);
  sprintf(out, "%.3f\t%.3f", *xd1, *xd2);
  for(int ia=0; ia<Na; ia++){
    Ave_se(Acc1+ia, Acc2+ia, num);    
    sprintf(tmp, "\t%.3f\t%.2f", Acc1[ia], Acc2[ia]);
    strcat(out, tmp);
    p= num_signi[ia]/(float)num;
    if(p > 1){
      printf("Error n_significant= %d  n=%d  alpha=%.2f\n",
	     num_signi[ia], num, alpha_val[ia]);
    }
    sprintf(tmp, "\t%.3f\t%.2f", p, sqrt(p*(1-p)/num));
    strcat(out, tmp);
  }
  Ave_se(t1, t2, num);
  p=*Sign1/(*Out1);
  sprintf(tmp, "\t%.3f\t%.2f\t%.3f\t%.2f",
	  *t1, *t2, p, sqrt(p*(1-p)/(*Out1)));
  strcat(out, tmp);
  p=*TV1/(*Out1+*TV1);
  sprintf(tmp, "\t%.3f\t%.2f\t%.1f\t%d",
	  p, sqrt(p*(1-p)/(*Out1+*TV1)), (*Out1)/num, num);
  strcat(out, tmp);
}

void Ave_se(double *x1, double *x2, int n)
{
  if(n==0)return;
  if(n==1){*x2=0; return;}
  *x1/=n;
  *x2=(*x2/n-(*x1)*(*x1))/(n-1);
  if(*x2<-0.00000001){
    printf("ERROR in st.dev., s^2= %f\n", *x2);
    return;
  }
  *x2=sqrt(*x2);
}

float Fit_alpha(float *thr,
		float **Div, float **CV, float **t, int **nout,
		float **Seq_Id, float SI_thr, int N,
		FILE *file_out)
{
  float coeff=1.1; // |CV|~c*d^alpha, Threshold=coeff*c
  float d_thr=0.3; // Fit made for D < d_thr <D>
  float a_max=1.0; // Maximum allowed value of alpha 
  float thr_def=0.15; // default value of threshold
  int i,j; double D1=0; int Np=0;
  for(i=0; i<N; i++){
    for(j=0; j<i; j++){
      if(Seq_Id[i][j]<SI_thr)continue;
      D1+=Div[i][j]; Np++;
    }
  }
  D1/=(float)Np;
  D1*=d_thr; 
  double X1=0, X2=0, Y1=0, Y2=0, XY=0, W1=0; int m=0;
  for(i=0; i<N; i++){
    for(j=0; j<i; j++){
      if((nout[i][j]<1)||(Div[i][j]<=0)||(Div[i][j]>D1)||(CV[i][j]==0))
	continue;
      if(Seq_Id[i][j]<SI_thr)continue;
      float x=log(Div[i][j]);
      float y=log(fabs(CV[i][j])); ///fabs(CV[i][j])
      //float w=t[i][j]; w*=w; // |<d(A,C)-d(B,C)>|/S.E.M.   
      //X1+=w*x; X2+=w*x*x; Y1+=w*y; Y2+=w*y*y; XY+=w*x*y; W1+=w;
      X1+=x; X2+=x*x; Y1+=y; Y2+=y*y; XY+=x*y; W1++;
      m++;
    }
  }
  double C11=X2;
  XY=XY*W1-X1*Y1;
  X2=X2*W1-X1*X1;
  Y2=Y2*W1-Y1*Y1;
  double alpha=XY/X2;
  double c=(Y1-alpha*X1)/W1;
  c=exp(c); *thr=coeff*c;
  if(file_out){
    double r=XY/sqrt(X2*Y2); // correlation coefficient
    // Error
    double Error=Y2*(1-r*r)/W1; Error/=(W1-2);
    double det=X2; // determinant of the covariance matrix
    double err_a=sqrt(Error*W1/det);    //(C^(-1))_11= W1/det
    double err_thr=sqrt(Error*C11/det); //(C^(-1))_00= sum_X^2/det
    err_thr*=c;
    fprintf(file_out, "Fit: |CV|~c*D^alpha D < %.2f<D>=%.2f %d pairs\n",
	    d_thr, D1, m);
    fprintf(file_out, "Correlation coefficient: %.3f\n", r);
    fprintf(file_out, "alpha= %.3f err= %.3f\n", alpha, err_a);
    fprintf(file_out, "coeff= %.3f err= %.3f\n", c, err_thr);
    fprintf(file_out, "threshold= %.3f (%.3f*c) err= %.3f\n",
	    *thr, coeff, coeff*err_thr);
    if(alpha>a_max){
      fprintf(file_out, "WARNING, setting alpha=%.1f max. allowed\n", a_max);
      fprintf(file_out, "WARNING, setting threshold= %.2f\n", thr_def);
    }
    fprintf(file_out, "\n");
  }
  if(alpha>a_max){
    alpha=a_max; *thr=thr_def;
  }
  return(alpha);
}

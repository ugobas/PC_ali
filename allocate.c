#include <stdlib.h>

float  **Allocate_mat2_f(int n1, int n2);
double **Allocate_mat2_d(int n1, int n2);

void Empty_matrix_f(float **matrix, int N){
  if(matrix==NULL)return;
  for(int i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}

void Empty_matrix_d(double **matrix, int N){
  if(matrix==NULL)return;
  for(int i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}


void Empty_matrix_i(int **matrix, int N){
  if(matrix==NULL)return;
  for(int i=0; i<N; i++)free(matrix[i]);
  free(matrix);
}

float **Allocate_mat2_f(int n1, int n2){
  float **matrix=malloc(n1*sizeof(float *));
  for(int i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(float));
    for(int j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}


int ** Allocate_mat2_i(int n1, int n2){
  int **matrix=malloc(n1*sizeof(int *));
  for(int i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(int));
    for(int j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}

double **Allocate_mat2_d(int n1, int n2){
  double **matrix=malloc(n1*sizeof(double *));
  for(int i=0; i<n1; i++){
    matrix[i]=malloc(n2*sizeof(double));
    for(int j=0; j<n2; j++)matrix[i][j]=(double)0;
  }
  return(matrix);
}

/*            Fortran             */

float **Allocate_mat2_f_fortran(int n1, int n2){
  int i, j;
  float **matrix=malloc(n1*sizeof(float *));
  matrix[0]=malloc((n1*n2)*sizeof(float));
  for(i=1; i<n1; i++){
    matrix[i]=matrix[i-1]+n2;
    for(j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}

int **Allocate_mat2_i_fortran(int n1, int n2){
  int **matrix=malloc(n1*sizeof(int *));
  matrix[0]=malloc((n1*n2)*sizeof(int));
  for(int i=1; i<n1; i++){
    matrix[i]=matrix[i-1]+n2;
    for(int j=0; j<n2; j++)matrix[i][j]=0;
  }
  return(matrix);
}

void Empty_matrix_f_fortran(int N, float **matrix){
  if(matrix==NULL)return;
  free(matrix[0]);
  free(matrix);
}

void Empty_matrix_i_fortran(int N, int **matrix){
  if(matrix==NULL)return;
  free(matrix[0]);
  free(matrix);
}

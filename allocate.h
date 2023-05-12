int ** Allocate_mat2_i(int n1, int n2);
float  **Allocate_mat2_f(int n1, int n2);
double **Allocate_mat2_d(int n1, int n2);
void Empty_matrix_i(int **matrix, int N);
void Empty_matrix_f(float **matrix, int N);
void Empty_matrix_d(double **matrix, int N);

float **Allocate_mat2_f_fortran(int n1, int n2);
int **Allocate_mat2_i_fortran(int n1, int n2);
void Empty_matrix_f_fortran(int N, float **matrix);
void Empty_matrix_i_fortran(int N, int **matrix);

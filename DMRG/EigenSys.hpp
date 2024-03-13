#include <complex>

using namespace std;

#define TRUE	1
#define FALSE	0

/////////////////////////////////////////////////////////////////////////////
// SparseMatrix Globals

// Minimal EigenValue and EigenVector
void      Davidson(int, SparseMatrix<double>&, Vector<double>&, Vector<double>&, int, double);

/////////////////////////////////////////////////////////////////////////////
// C subroutines

void tred2_c(double **a, int n, double d[], double e[]);

void tqli_c(double d[], double e[], int n, double **z);

void sort_c(int n, int SMALL, double d[], double **v);

void sort_qn_c(int n, int SMALL, double d[], double *v, double qn[], long index[], int index2[]);

void sort_qn_2_c(int n, int SMALL, double d[], double *v, int qn_up[], int qn_down[], long index[], int index2[]);

void sort_qn_3_c(int n, int SMALL, double d[], double *v, int **qn_LA, double qn_d[], long index[], int index2[]);
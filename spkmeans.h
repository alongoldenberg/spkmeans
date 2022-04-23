/**
* Created by Alon Goldenberg & Ofir Nissan.
*/

#ifndef SPKMEANS_SPKMEANS_H
#define SPKMEANS_SPKMEANS_H

#define PY_SSIZE_T_CLEAN

#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define MAX_ROTATIONS 100
#define MAX_ITER 300
#define EPSILON pow(10,-5)

typedef struct {
    /* A structure containing index and value */
    int i;
    double v;
} indexed_double;

double **allocate_data(int n, int d);
double **parse_file(char *filename, int *n, int *d);
double **weight_adj_matrix(double **datapoints, int n, int d);
double *diagonal_degree_matrix(double **weights, char sq, int n);
double **normalized_laplacian(double **weights, const double *degree, int n);
double rotate_jacobian(double **a, double **a_tag, double **p, int n);
double **jacobi_function(double **a, double eps, int n);
double calc_c(double t);
double calc_t(double **m, int i, int j);
double sign(double x);
void max_off_diagonal(double **m, int n, int *max_i, int *max_j);
void calc_a_tag(double **a, double **a_tag,
                double s, double c, int i, int j, int n);
void calc_p(double **p, int i, int j, double s, double c, int n);
void copy_matrix(double **m, double **c, int n);
double **identity_matrix(int k);
double *get_diagonal(double **m, int n);
double **degree_to_diagonal_matrix(const double *degree, int n);
int eigengap_hueuristic(const double *eigenvaleus, int n);
double **kmeans(double **datapoints, double **centroids,
                int k, int max_iter, double epsilon, int n, int d);
double update_centroids(double **centroids,
                        double **datapoints, int , int n, int d);
void update_cumulative_sums(const double *arr, double *cumulative_sum, int d);
void sort_eigenvalues_and_vectors(const double *eigenvalues,
                                  double **eigenvectors, double * s_eigenvalues,
                                  double ** s_eigenvectors, int n);
double **calculate_T(double **eigenvectors, int k, int n);
double **spectral_clustrering(double **datapoints, int n, int d, int k);
void transpose(double** A, int n);

#endif

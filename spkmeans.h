/**
* Created by Alon Goldenberg & Ofir Nissan.
*/

#ifndef SPKMEANS_SPKMEANS_H
#define SPKMEANS_SPKMEANS_H

#define PY_SSIZE_T_CLEAN

#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define MAX_ROTATIONS 100
#define MAX_ITTER 300
#define EPSILON exp(-15)


static double **allocate_data(int n, int d);
double **parse_file(char *filename, int *n, int *d);
static double **weight_adj_matrix(double **datapoints, int n, int d);
static double *diagonal_degree_matrix(double **weights, char sq, int n);
static double **normalized_laplacian(double **weights, const double *degree, int n);
static double rotate_jacobian(double **a, double **a_tag, double **p, int n);
static double **jacobi_function(double **a, double eps, int n);
static double calc_c(double t);
static double calc_t(double **m, int i, int j);
static double sign(double x);
static int *max_off_diagonal(double **m, int n);
static void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j, int n);
static void calc_p(double **p, int i, int j, double s, double c, int n);
static void copy_matrix(double **m, double **c, int n);
static double **identity_matrix(int k);
static double *get_diagonal(double **m, int n);
static double **degree_to_diagonal_matrix(const double *degree, int n);
static int eigengap_hueuristic(const double *eigenvaleus, int n);
static double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon, int n, int d);
static double update_centroids(double **centroids, double **datapoints, int , int n, int d);
static void update_cumulative_sums(const double *arr, double *cumulative_sum, int d);
static void sort_eigenvalues_and_vectors(const double *eigenvalues, double **eigenvectors,
                                         double * s_eigenvalues, double ** s_eigenvectors, int n);
static double **calculate_T(double **eigenvectors, int k, int n);
static double **calculate_k(double **datapoints, int n, int d);
static double **spectral_clustrering(double **datapoints, int n, int d);



#endif

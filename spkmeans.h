

#ifndef SPKMEANS_SPKMEANS_H
#define SPKMEANS_SPKMEANS_H

#define PY_SSIZE_T_CLEAN

#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define MAX_ROTATIONS 100
#define MAX_ITTER 300
#define EPSILON exp(-15)


double **allocate_data(int n, int d);
double **parse_file(char *filename);
double **parse_file(char *filename);
double **weight_adj_matrix(double **datapoints);
double *diagonal_degree_matrix(double **weights, char sq);
double **normalized_laplacian(double **weights, const double *degree);
double rotate_jacobian(double **a, double **a_tag, double **p);
double **jacobi_function(double **a, double eps);
double calc_c(double t);
double calc_t(double **m, int i, int j);
double sign(double x);
int *max_off_diagonal(double **m);
void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j);
void calc_p(double **p, int i, int j, double s, double c);
void copy_matrix(double **m, double **c);
double **identity_matrix(int k);
double *get_diagonal(double **m);
double **degree_to_diagonal_matrix(double *degree);
int eigengap_hueuristic(const double *eigenvaleus);
int calculate_k(double **datapoints);
double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon);
double update_centroids(double **centroids, double **datapoints, int k);
void update_cumulative_sums(double *arr, double *cumulative_sum);
void sort_eigenvalues_and_vectors(const double *eigenvalues, double **eigenvectors,
                                         double * s_eigenvalues, double ** s_eigenvectors);
double **calculate_T(double **eigenvectors, int k);


#endif 

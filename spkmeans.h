//
// Created by Alon Goldenberg on 20/02/2022.
//

#ifndef SPKMEANS_SPKMEANS_H
#define SPKMEANS_SPKMEANS_H

#define PY_SSIZE_T_CLEAN
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define MAX_ROTATIONS 100
#define MAX_ITTER 300
#define EPSILON exp(-15)


static double **allocate_data(int n, int d);
double **parse_file(char *filename);
double **parse_file(char *filename);
static double **weight_adj_matrix(double **datapoints);
static double *diagonal_degree_matrix(double **weights, char sq);
static double **normalized_laplacian(double **weights, const double *degree);
static double rotate_jacobian(double **a, double **a_tag, double **p);
static double **jacobi_function(double **a, double eps);
static double calc_c(double t);
static double calc_t(double **m, int i, int j);
static double sign(double x);
static int *max_off_diagonal(double **m);
static void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j);
static void calc_p(double **p, int i, int j, double s, double c);
static void copy_matrix(double **m, double **c);
static double **identity_matrix(int k);
static double *get_diagonal(double **m);
static double **degree_to_diagonal_matrix(double *degree);
static int eigengap_hueuristic(const double *eigenvaleus);
static int calculate_k(double **datapoints);
static double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon);
static double update_centroids(double **centroids, double **datapoints, int k);
static void update_cumulative_sums(double *arr, double *cumulative_sum);
static void sort_eigenvalues_and_vectors(const double *eigenvalues, double **eigenvectors,
                                         double * s_eigenvalues, double ** s_eigenvectors);
static double **calculate_T(double **eigenvectors, int k);


#endif //SPKMEANS_SPKMEANS_H

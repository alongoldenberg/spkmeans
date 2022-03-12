//
// Created by Alon Goldenberg on 20/02/2022.
//

#ifndef SPKMEANS_SPKMEANS_H
#define SPKMEANS_SPKMEANS_H

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define MAX_ROTATIONS 100

static double **allocate_data(int n, int d);
static double **weight_adj_matrix(double **datapoints);
static double *diagonal_degree_matrix(double **weights);
static double **normalized_laplacian(double **weights, double *degree);
static double sign(double x);
static double calc_c(double t);
static double calc_t(double **m, int i, int j);
static double **copy_matrix(double **m, double **c);
static int *max_off_diagonal(double **m);
static double calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j);

#endif //SPKMEANS_SPKMEANS_H

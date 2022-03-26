//
// Created by Alon Goldenberg & Ofir Nissan.
//

#include "spkmeans.h"
#include "utils.c"

int n, d;

// TODO: Make sure to free all matrix

static double **weight_adj_matrix(double **datapoints){
/**
 * Calculate Weights Matrix
 *
 * @param datapoints - n*d matrix
 * @return Weights matrix - n*n matrix where W[i,j] = weight between datapoint i and j
 */
    int i, j; double d, w; double ** weights;
    weights = allocate_data(n, n);
    for (i=0;i<n;i++){
        for(j=0;j<=i;j++){
            if(i==j){
                weights[i][i] = 0;
                continue;
            }
            d = distance(datapoints[i], datapoints[j], d) / 2;
            w = exp(-d);
            weights[i][j] = w;
            weights[j][i] = w;
        }
    }
    return weights;
}


static double *diagonal_degree_matrix(double **weights){
/**
 * Calculate degree matrix
 *
 * @param weights - n*n weights matrix
 * @return degree array - n sized array where d[i] = degree of datapoint i
 */
    // Note: this function returns a 1D array but is considered as 2D with A[i][i] = m[i] with 0 else.
    double sum; double *diagonal; int i,j;
    diagonal = (double*) calloc(n, sizeof(double));
    for (i=0;i<n;i++){
        sum = 0;
        for (j=0;j<n;j++){
            sum+=weights[i][j];
        }
        sum = 1 / (sqrt(sum));
        diagonal[i] = sum;
    }
    return diagonal;
}

static double **normalized_laplacian(double **weights, double *degree){
/**
 * Calculate normalized laplacian matrix
 *
 * @param weights - n*n weights matrix
 * @param degree - n size degree array
 * @return normlized laplacian matrix - n*n matrix represnting the normalized laplacina
 */
    int i,j; double **lap; double x;
    lap = allocate_data(n, n);
    for (i=0;i<n;i++){
        for(j=0;j<=i;j++){
            if(i==j){
                lap[i][i] = 1;
                continue;
            }
            x = -(weights[i][j]*degree[i]*degree[j]);
            lap[i][j] = x;
            lap[j][i] = x;
        }
    }
    return lap;
}

static double *jacobi(double **m, int eps){
    int i,j, rotations=0; double *max;
    max = max_off_diagonal(m);
    i = max[0], j=max[1];

}

static double calc_c(double t){
// calculates value of c
    double c;
    c =  1 / sqrt((pow(t, 2) + 1));
    return c;
}

static double calc_t(double **m, int i, int j){
    // calculates value of t
    double theta, t;
    theta = (m[j][j] - m[i][i]) / m[i][j];
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    return t;
}

static double sign(double x) {
    // given double returns it's sign. Note: sign(0) = 1.
    double s;
    if(x == 0){
        return 1;
    }
    s = (double)(x > 0) - (x < 0);
    return s;
}

static int *max_off_diagonal(double **m){
/**
 * Finds the indices of the maximum double in a matrix
 *
 * @param m - n*n matrix
 * @return [i,j] - 2-array of int represting row and column of the maximal double
 */
    double *res; double max = 0; int i,j;
    res =  (int*) calloc(2, sizeof(int));
    for(i=0;i<n;i++){
        for(j=0;j<i;j++){
            if (m[i][j] > max){
                max = m[i][j];
                res[0] = i;
                res[1] = j;
            }
        }
    }
    return res;
}

static double calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j){
/**
 * Calculates A' and returns the off-diagonal difference between A and A'
 *
 * @param a - n*n matrix representing a
 * @param a_tag - n*n empty matrix to be filled with values
 * @param s, c - values for the transformation
 * @param i,j - indices of the maximal number in A
 * @return off - off-diagonal difference between A and A'
 */    int r; double diff=0;
    copy_matrix(a, a_tag);
    for(r=0; r<n; r++){
        if (r!=i && r!=j){
            a_tag[r][i] = c*a[r][i] - s*a[r][j];
            a_tag[i][r] = c*a[r][i] - s*a[r][j];
            a_tag[r][j] = c*a[r][j] + s*a[r][i];
            a_tag[j][r] = c*a[r][j] + s*a[r][i];
            diff += 2 * (pow(a[r][i], 2) - pow(a_tag[r][i], 2));
            diff += 2 * (pow(a[r][j], 2) - pow(a_tag[r][j], 2));
        }
    }
    a_tag[i][i] = pow(c,2) * a[i][i] + pow(s,2) * a[j][j] - 2*s*c*a[i][j];
    a_tag[j][j] = pow(s,2) * a[i][i] + pow(c,2) * a[j][j] + 2*s*c*a[i][j];
    a_tag[i][j] = 0;
    return diff;
}

static double **copy_matrix(double **m, double **c) {
    // Deep copy a n*n matrix
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = m[i][j];
        }
    }
}


static double **degree_to_diagonal_matrix(double *degree){
    double **diagonal_degree_matrix_res;
    diagonal_degree_matrix_res = allocate_data(n, n);
    int i,j;
    for (i=0 ; i<n; i++){
        for (j=0 ; j<n; j++){
            diagonal_degree_matrix_res[i][j] = degree[i];
        }
    }
    return diagonal_degree_matrix_res;
}


int main(int argc, char *argv[]){
    enum Goal {wam, ddg, lnorm, jacobi} goal;
    char *file_name;
    double **datapoints;
    int k;
    double **weight_adj_matrix_res; double *degree; double **lap_res;
    double **diagonal_degree_matrix_res;
    if (argc == 4) {
        k = atoi(argv[1]);
        goal = argv[2];
        file_name = argv[3];
    }
    else{
        printf("Invalid Input!");
        return 1;
    }
    datapoints = parse_file(file_name);
    switch (goal)
    {
        case wam:
            weight_adj_matrix_res = weight_adj_matrix(datapoints);
            print_data(weight_adj_matrix_res);
            break;
        case ddg:
            degree = diagonal_degree_matrix(weight_adj_matrix_res);
            diagonal_degree_matrix_res = degree_to_diagonal_matrix(degree);
            print_data(diagonal_degree_matrix_res);
            break;
        case lnorm:
            lap_res = normalized_laplacian(weight_adj_matrix_res, degree);
            print_data(lap_res);
            break;
        case jacobi:
        //need to complete;
        default:
            printf("Invalid Input!");
        return 1;
    }  
    return 0;
}

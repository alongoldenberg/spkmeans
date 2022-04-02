//
// Created by Alon Goldenberg & Ofir Nissan.
//

#include "spkmeans.h"
#include "utils.c"

int n, d;

// TODO: Make sure to free all matrix
// TODO: Check that all calloc succeeded.

static double **weight_adj_matrix(double **datapoints){
/**
 * Calculate Weights Matrix
 *
 * @param datapoints - n*d matrix
 * @return Weights matrix - n*n matrix where W[i,j] = weight between datapoint i and j
 */
    int i, j; double dist, w; double ** weights;
    weights = allocate_data(n, n);
    for (i=0;i<n;i++){
        for(j=0;j<=i;j++){
            if(i==j){
                weights[i][i] = 0;
                continue;
            }
            dist = distance(datapoints[i], datapoints[j], d) / 2;
            w = exp(-dist);
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
 * @return normalized laplacian matrix - n*n matrix representing the normalized laplacian
 */
    int i,j; double **lap; double x;
    lap = allocate_data(n, n);
    for (i=0;i<n;i++){
        for(j=0;j<=i;j++){
            if(i==j){
                lap[i][i] = 1;
                continue;
            }
            x = -(weights[i][j] * degree[i] * degree[j]);
            lap[i][j] = x;
            lap[j][i] = x;
        }
    }
    return lap;
}

static double rotate_jacobian(double **a, double **a_tag, double **p) {
    int i, j; double s, c, t; int *max; double max_value;
    max = max_off_diagonal(a);
    if (max==0){
        // A is already diagonal
        return 0;
    }
    i = max[0], j=max[1];
    max_value = a[i][j];
    t = calc_t(a, i, j);
    c = calc_c(t);
    s = c * t;
    calc_a_tag(a, a_tag, s, c, i, j);
    calc_p(p, i, j, s, c);
    copy_matrix(a_tag, a);
    return 2*(pow(max_value, 2));
}

static double **jacobi(double **a){
/**
 * Preform the Jacobian algorithm.
 *
 * @param a - n*n real symmetric metrix
 * @param eps - tolerance for the difference between rotation matrices
 * @return **vectors  -  2D array with the eighenvectors as columns
 */
    int k; double off_diagonal;
    double **a_tag, **p;
    a_tag = allocate_data(n, n);
    copy_matrix(a, a_tag);
    p = identity_matrix(n);
    for(k=0;k<MAX_ROTATIONS;k++){
        off_diagonal = rotate_jacobian(a, a_tag, p);
        if(off_diagonal <= EPSILON) {
            break;
        }
    }
  return p;
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
 * @return [i,j] - 2-array of int representing row and column of the maximal double
 */
    int *res; double max = 0; int i,j;
    res =  (int*) calloc(2, sizeof(int));
    for(i=0;i<n;i++){
        for(j=0;j<i;j++){
            if (fabs(m[i][j]) > fabs(max)){
                max = m[i][j];
                res[0] = i;
                res[1] = j;
            }
        }
    }
    return res;
}

static void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j){
/**
 * Calculates A' and returns the off-diagonal difference between A and A'
 *
 * @param a - n*n matrix representing a
 * @param a_tag - n*n empty matrix to be filled with values
 * @param s, c - values for the transformation
 * @param i,j - indices of the maximal number in A
 * @return off - off-diagonal difference between A and A', based on
 */
    int r;
    copy_matrix(a, a_tag);
    for(r=0; r<n; r++){
        if (r!=i && r!=j){
            a_tag[r][i] = c*a[r][i] - s*a[r][j];
            a_tag[i][r] = c*a[r][i] - s*a[r][j];
            a_tag[r][j] = c*a[r][j] + s*a[r][i];
            a_tag[j][r] = c*a[r][j] + s*a[r][i];
        }
    }
    a_tag[i][i] = pow(c,2) * a[i][i] + pow(s,2) * a[j][j] - 2*s*c*a[i][j];
    a_tag[j][j] = pow(s,2) * a[i][i] + pow(c,2) * a[j][j] + 2*s*c*a[i][j];
    a_tag[i][j] = 0;
}

static void calc_p(double **p, int i, int j, double s, double c){
/**
 * Calculates Rotation Matrices in place using
 *
 * @param p - n*n matrix representing current rotation matrix
 * @param s, c - values for the transformation
 * @param i,j - indices of the maximal number in A
 */
    double tau, temp; int r;
    tau = s / (1.0 + c);
    for(r=0; r<n; r++){
        temp = p[r][i];
        p[r][i] = temp - s*(p[r][j] + tau*p[r][i]);
        p[r][j] = p[r][i] + s*(temp - tau*p[r][j]);
    }
}

static void copy_matrix(double **m, double **c) {
    // Deep copy a n*n matrix
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = m[i][j];
        }
    }
}

static double **identity_matrix(int k) {
    // return an N*N identity matrix
    double **identity; int i;
    identity = allocate_data(n, n);
    for (i = 0; i < k; i++) {
        identity[i][i] = 1;
    }
    return identity;
}

static double *get_diagonal(double **m) {
    // return an N*N identity matrix
    double *diagonal; int i;
    diagonal = (double *)calloc(n, sizeof(double ));
    for (i = 0; i < n; i++) {
        diagonal[i] = m[i][i];
    }
    return diagonal;
}



static double **degree_to_diagonal_matrix(double *degree){
    // Use to convert the 1D array to 2D array
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

int compare( const void* a, const void* b)
{
    // simple compare function.
    int int_a = * ( (int*) a );
    int int_b = * ( (int*) b );

    if ( int_a == int_b ) {
        return 0;
    }
    else if ( int_a < int_b )
    {return -1;}
    return 1;
}


static int eigengap_hueuristic(double *eigenvaleus){
/**
 * Preform the eigengap heuristic: given a set of eigenvalues return max idx difference of sorted eigenvalues
 *
 * @param eigenvalues - sorted eigenvalues array
  *@result max_diff  - the maximal difference idx in 1-(n/2) smallest eigenvalues
 */
    int max_diff_idx, i;
    double max_diff = 0, diff;
    for (i = 0; i < n / 2; i++) {
        diff = eigenvaleus[i + 1] - eigenvaleus[i];
        if (diff > max_diff) {
            max_diff = diff;
            max_diff_idx = i;
        }
    }
    return max_diff_idx + 1;
}

static int calculate_k(double **datapoints){
/**
 * Preform the full spectral k-means algorithm, return optimal k.
 *
 * @param datapoints - 2D array of datapoints size n*d.
 * @invariant - n and d are initialized.
  *@result k - optimal k for k-means algorithm.
 */    
    double **weights, *degree, **laplacian, *eigenvalues;
    int k;
    weights = weight_adj_matrix(datapoints);
    degree = diagonal_degree_matrix(weights);
    laplacian = normalized_laplacian(weights, degree);
    jacobi(laplacian); // No need to return as laplacian is diagonalized in-place.
    eigenvalues = get_diagonal(laplacian);
    qsort(eigenvalues, n, sizeof (double), compare);
    k = eigengap_hueuristic(eigenvalues);
    return k;
}
/**
 For running as main:
 * */


//int main(int argc, char *argv[]){
//    double **w;
//    w = weight_adj_matrix(datapoints);
//    print_matrix(w, n, d);
//
//
//    return 0;
//}

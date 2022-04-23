/**
* Created by Alon Goldenberg & Ofir Nissan.
*/

#include "spkmeans.h"
#include "utils.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

void transpose(double** A, int n){
    int i; int j; double tmp;
    for (i = 0; i < n - 1; i++){
        for (j = i + 1; j < n; j++){
            tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    }
}

double **parse_file(char *filename, int* n, int* d){
    /* Allocate memory and put data from file to a 2d array */
    FILE *data = NULL; double **arr; double num; int i,j;

    data = fopen(filename, "r");
    if(data == NULL){
        print_invalid_input();
    }
    *n = count_lines(data);
    *d = count_d(data);
    arr = allocate_data(*n, *d);
    rewind(data);
    for(i = 0; i<*n; i++){
        for(j=0; j<*d; j++) {
            fscanf(data, "%lf,", &num);
            arr[i][j] = num;
        }
    }
    fclose(data);
    return arr;
}



double **weight_adj_matrix(double **datapoints, int n, int d){
/**
 * Calculate Weights Matrix
 *
 * @param datapoints - n*d matrix
 * @return Weights matrix - n*n matrix where W[i,j] = weight between datapoint i and j
 */
    int i, j; double dist, w; double ** weights;
    weights = allocate_data(n, n);
    for (i=0;i<n;i++){
        for(j=0;j<i;j++){
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



double *diagonal_degree_matrix(double **weights, char sq, int n){
/**
 * Calculate degree matrix
 *
 * @param weights - n*n weights matrix
 * @param sq - if !=0, return D**-0.5
 * @return degree array - n sized array where d[i] = degree of datapoint i
 * Note: this function returns a 1D array but is considered as 2D with A[i][i] = m[i] with 0 else.
 */
    double sum; double *diagonal; int i,j;
    diagonal = (double*) calloc(n, sizeof(double));
    if (diagonal == NULL) {
        print_error();
    }
    for (i=0;i<n;i++){
        sum = 0;
        for (j=0;j<n;j++){
            sum+=weights[i][j];
        }
        if(sq) {
            sum = 1 / (sqrt(sum));
        }
        diagonal[i] = sum;
    }
    return diagonal;
}


double **normalized_laplacian(double **weights, const double *degree, int n){
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


double rotate_jacobian(double **a, double **a_tag, double **p, int n) {
    int i, j; double s, c, t; double max_value;
    i=0;
    j=1;
    max_off_diagonal(a, n, &i, &j);
    if (a[i][j]==0){
        return 0;
    }
    max_value = a[i][j];
    t = calc_t(a, i, j);
    c = calc_c(t);
    s = c * t;
    calc_a_tag(a, a_tag, s, c, i, j, n);
    calc_p(p, i, j, s, c, n);
    copy_matrix(a_tag, a, n);
    return 2*(pow(max_value, 2));
}


double **jacobi_function(double **a, double eps, int n){
/**
 * Preform the Jacobian algorithm.
 *
 * @param a - n*n real symmetric matrix
 * @param eps - tolerance for the difference between rotation matrices
 * @return **vectors  -  2D array with the eighenvectors as columns
 */
    int k; double off_diagonal;
    double **a_tag, **p;
    a_tag = allocate_data(n, n);
    copy_matrix(a, a_tag, n);
    p = identity_matrix(n);
    for(k=0;k<MAX_ROTATIONS;k++){
        off_diagonal = rotate_jacobian(a, a_tag, p, n);
        if(off_diagonal <= eps) {
            break;
        }
    }
    free(a_tag[0]);
    free(a_tag);
    return p;
}


double calc_c(double t){
/**calculates value of c*/
    double c;
    c =  1 / sqrt((pow(t, 2) + 1));
    return c;
}


double calc_t(double **m, int i, int j){
/** calculates value of t*/
    double theta, t;
    theta = (m[j][j] - m[i][i]) / (2*m[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    return t;
}


double sign(double x) {
/** given double returns it's sign. Note: sign(0) = 1. */
    double s;
    if(x == 0){
        return 1;
    }
    s = (double)(x > 0) - (x < 0);
    return s;
}

void max_off_diagonal(double **m, int n, int *max_i, int *max_j){
/**
 * Finds the indices of the maximum double in a matrix
 *
 * @param m - n*n matrix
 */
    double max = 0; int i,j;
    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            if (fabs(m[i][j]) > fabs(max)){
                max = m[i][j];
                *max_i = i;
                *max_j = j;
            }
        }
    }
}


void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j, int n){
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
    copy_matrix(a, a_tag, n);
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
    a_tag[j][i] = 0;
}


void calc_p(double **p, int i, int j, double s, double c, int n){
/**
 * Calculates Rotation Matrices in place using heuristics
 *
 * @param p - n*n matrix representing current rotation matrix
 * @param s, c - values for the transformation
 * @param i,j - indices of the maximal number in A
 */
    double temp; int r;
    for(r=0; r<n; r++){
        temp = p[r][i] * s + p[r][j] * c;
        p[r][i] = p[r][i] * c - p[r][j] * s;
        p[r][j] = temp;
    }
}


void copy_matrix(double **m, double **c, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = m[i][j];
        }
    }
}


double **identity_matrix(int n) {
    double **identity; int i;
    identity = allocate_data(n, n);
    for (i = 0; i < n; i++) {
        identity[i][i] = 1;
    }
    return identity;
}


double *get_diagonal(double **m, int n) {
    double *diagonal; int i;
    diagonal = (double *)calloc(n, sizeof(double ));
    if (diagonal == NULL) {
        print_error();
    }
    for (i = 0; i < n; i++) {
        diagonal[i] = m[i][i];
    }
    return diagonal;
}


double **degree_to_diagonal_matrix(const double *degree, int n){
    double **diagonal_degree_matrix_res;
    int i;
    diagonal_degree_matrix_res = allocate_data(n, n);
    for (i=0 ; i<n; i++){
            diagonal_degree_matrix_res[i][i] = degree[i];
    }
    return diagonal_degree_matrix_res;
}


int idx_cmp( const void* a, const void* b)
{

    indexed_double *x = (indexed_double *) a;
    indexed_double *y = (indexed_double *) b;

    if ( (*x).v == (*y).v ) {
        return 0;
    }
    else if ( (*x).v < (*y).v )
    {return -1;}
    return 1;
}



int eigengap_hueuristic(const double *eigenvaleus, int n){
/**
 * Preform the eigengap heuristic: given a set of eigenvalues return max idx difference of sorted eigenvalues
 *
 * @param eigenvalues - sorted eigenvalues array
  *@result max_diff  - the maximal difference idx in 1-(n/2) smallest eigenvalues
 */
    int max_diff_idx = 0, i;
    double max_diff = 0, diff;
    for (i = 0; i <= n / 2; i++) {
        diff = fabs(eigenvaleus[i + 1] - eigenvaleus[i]);
        if (diff > max_diff) {
            max_diff = diff;
            max_diff_idx = i;
        }
    }
    return max_diff_idx + 1;
}

void sort_eigenvalues_and_vectors(const double *eigenvalues, double **eigenvectors,
                                         double *s_eigenvalues, double **s_eigenvectors, int n){
    indexed_double *eigenvalues_idx;
    int i;
    eigenvalues_idx = (indexed_double *) calloc(n, sizeof (indexed_double));
    if (eigenvalues_idx == NULL) {
        print_error();
    }
    for(i = 0; i<n; i++){
        eigenvalues_idx[i].i = i;
        eigenvalues_idx[i].v = eigenvalues[i];
    }
    qsort(eigenvalues_idx, n, sizeof (eigenvalues_idx[0]), idx_cmp);
    transpose(eigenvectors, n);
    for(i = 0; i<n; i++){
        s_eigenvalues[i] = eigenvalues[eigenvalues_idx[i].i];
        s_eigenvectors[i] = eigenvectors[eigenvalues_idx[i].i];
    }
    transpose(s_eigenvectors, n);
    free(eigenvalues_idx);
}


double **calculate_T(double **eigenvectors, int k, int n){
    double **T, *norms, sum;
    int i,j;
    T = allocate_data(n, k);
    norms = (double *) calloc(n, sizeof (double));
    if (norms == NULL) {
        print_error();
    }
    for (i=0;i<n;i++){
        sum = 0;
        for(j=0;j<k;j++){
            sum+= pow(eigenvectors[i][j], 2);
        }
        norms[i] = sqrt(sum);
    }
    for (i=0;i<n;i++){
        for(j=0;j<k;j++){
            if(norms[i] == 0){
                T[i][j] = 0;
            }
            else {
                T[i][j] = (eigenvectors[i][j] / norms[i]);
            }
        }
    }
    free(norms);
    return T;
}



double **spectral_clustrering(double **datapoints, int n, int d, int k){
/**
 * Preform the full spectral k-means algorithm, return k first eigenvectors.
 *
 * @param datapoints - 2D array of datapoints size n*d.
 * @invariant - n and d are initialized.
  *@result T - requested matrix 
 */    
    double **weights, *degree, **laplacian, *eigenvalues, **eigenvectors,
            *s_eigenvalues, **s_eigenvectors, **T;
    s_eigenvalues = (double *) calloc(n, sizeof (double));
    if (s_eigenvalues == NULL) {
        print_error();
    }
    s_eigenvectors = allocate_data(n, n);
    if (s_eigenvectors == NULL) {
        print_error();
    }
    weights = weight_adj_matrix(datapoints, n, d);
    degree = diagonal_degree_matrix(weights, 1, n);
    laplacian = normalized_laplacian(weights, degree, n);

    eigenvectors = jacobi_function(laplacian, EPSILON, n);
    eigenvalues = get_diagonal(laplacian, n);
    sort_eigenvalues_and_vectors(eigenvalues, eigenvectors, s_eigenvalues, s_eigenvectors, n);
    if(k==0) {
        k = eigengap_hueuristic(s_eigenvalues, n);
    }
    T = calculate_T(s_eigenvectors, k, n);
    free(degree);
    free(weights[0]);
    free(weights);   
    free(laplacian[0]);
    free(laplacian);
    free(s_eigenvalues);   
    free(eigenvalues);
    free(eigenvectors[0]);
    free(eigenvectors);
    free(s_eigenvectors);
    return T;
}


/*kmeans from first and second exc:*/
double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon, int n, int d) {
    double max_change; int iterations=0;
    do{
        max_change = update_centroids(centroids, datapoints, k, n, d);
        iterations++;
    } while ((max_change >= epsilon) && (iterations < max_iter));
    return centroids;
}


double update_centroids(double **centroids, double **datapoints, int k, int n, int d){
    double **cumulative_sums; double *counters; int chosen_m_idx; double max_change=0; double *old_centroid; int i,j;
    counters = (double*) calloc(k*sizeof(double), sizeof(double));
    if (counters == NULL) {
        print_error();
    }
    cumulative_sums = allocate_data(k,d);
    for (i=0; i<n; i++){
        chosen_m_idx = 0;
        for (j=0; j<k; j++){
            if (distance(centroids[j], datapoints[i], d) < distance(centroids[chosen_m_idx], datapoints[i], d)){
                chosen_m_idx = j;
            }
        }
        update_cumulative_sums(datapoints[i], cumulative_sums[chosen_m_idx], d);
        counters[chosen_m_idx] += 1;
    }
    /* calculate new k centroids
     and calculate the maximum euclidean norm ||delta_mu||:*/
    old_centroid = (double*)calloc(n*d, sizeof(double));
    if (old_centroid == NULL){
        print_error();
    }
    for (i=0; i<k; i++){
        for (j=0; j<d; j++){
            old_centroid[j] = centroids[i][j];
            centroids[i][j] = cumulative_sums[i][j] / counters[i];
        }
        if ((distance(old_centroid, centroids[i], d)) > max_change){
            max_change = distance(old_centroid, centroids[i], d);
        }
    }
    free(counters);
    free(cumulative_sums[0]);
    free(cumulative_sums);
    free(old_centroid);
    return max_change;
}


void update_cumulative_sums(const double *arr, double *cumulative_sum, int d){
    int i;
    for (i=0; i<d; i++){
        cumulative_sum[i] += arr[i];
    }
}


int main(int argc, char *argv[]) {
    char *file_name, *goal;
    int n, d;
    double **eigen_vectors, **diagonal_degree_matrix_res, **lap_res, *degree,
        **datapoints, **weight_adj_matrix_res, *eigenvalues;
    if (argc != 3) {
        print_invalid_input();

    }
    goal = argv[1];
    file_name = argv[2];
    datapoints = parse_file(file_name, &n, &d);
    
    if(strcmp(goal, "wam") == 0) {
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        print_matrix(weight_adj_matrix_res, n, n);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
    }
    else if(strcmp(goal, "ddg") == 0) {

        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 0, n);
        diagonal_degree_matrix_res = degree_to_diagonal_matrix(degree, n);
        print_matrix(diagonal_degree_matrix_res, n, n);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
        free(diagonal_degree_matrix_res[0]);
        free(diagonal_degree_matrix_res);
        free(degree);
    }
    else if(strcmp(goal, "lnorm") == 0) {
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 1, n);
        lap_res = normalized_laplacian(weight_adj_matrix_res, degree, n);
        print_matrix(lap_res, n, n);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
        free(degree);
        free(lap_res[0]);
        free(lap_res);
    }
    else if(strcmp(goal, "jacobi") == 0) {
        if (n!=d){
            printf("Invalid Input\n");
            free(datapoints[0]);
            free(datapoints);
            return 1;
        }
        eigen_vectors = jacobi_function(datapoints, EPSILON, n);
        eigenvalues = get_diagonal(datapoints, n);
        print_eigenvalues(eigenvalues, n);
        print_matrix(eigen_vectors, n, n);
        free(eigenvalues);
        free(eigen_vectors[0]);
        free(eigen_vectors);
    }
    else{
        printf("Invalid Input\n");
        free(datapoints[0]);
        free(datapoints);
        return 1;
    }
    free(datapoints[0]);
    free(datapoints);
    return 0;
}


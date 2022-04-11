
#include "spkmeans.h"
#include "utils.c"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

int n, d;

typedef struct {
    int i;
    double v;
} indexed_double;



double **parse_file(char *filename){
    /* Allocate memory and put data from file to a 2d array */
    FILE *data = NULL; double **arr; double num; int i,j;

    data = fopen(filename, "r");
    if(data == NULL){
        print_invalid_input();
    }
    n = count_lines(data);
    d = count_d(data);
    arr = allocate_data(n, d);
    rewind(data);
    for(i = 0; i<n; i++){
        for(j=0; j<d; j++) {
            fscanf(data, "%lf,", &num);
            arr[i][j] = num;
        }
    }
    fclose(data);
    return arr;
}


double **weight_adj_matrix(double **datapoints){
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


double *diagonal_degree_matrix(double **weights, char sq){
/**
 * Calculate degree matrix
 *
 * @param weights - n*n weights matrix
 * @param sq - if !=0, returm D**-0.5
 * @return degree array - n sized array where d[i] = degree of datapoint i
 */
    /* Note: this function returns a 1D array but is considered as 2D with A[i][i] = m[i] with 0 else.*/
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

double **normalized_laplacian(double **weights, const double *degree){
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

double rotate_jacobian(double **a, double **a_tag, double **p) {
    int i, j; double s, c, t; int *max; double max_value;
    max = max_off_diagonal(a);
    i = max[0], j=max[1];
    if (a[i][j]==0){
        /* A is already diagonal*/
        return 0;
    }
    max_value = a[i][j];
    t = calc_t(a, i, j);
    c = calc_c(t);
    s = c * t;
    calc_a_tag(a, a_tag, s, c, i, j);
    calc_p(p, i, j, s, c);
    copy_matrix(a_tag, a);
    return 2*(pow(max_value, 2));
}


double **jacobi_function(double **a, double eps){
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
    copy_matrix(a, a_tag);
    p = identity_matrix(n);
    for(k=0;k<MAX_ROTATIONS;k++){
        off_diagonal = rotate_jacobian(a, a_tag, p);
        if(off_diagonal <= eps) {
            break;
        }
    }
    free(a_tag);
    return p;
}

double calc_c(double t){
/* calculates value of c*/
    double c;
    c =  1 / sqrt((pow(t, 2) + 1));
    return c;
}

double calc_t(double **m, int i, int j){
    /* calculates value of t*/
    double theta, t;
    theta = (m[j][j] - m[i][i]) / (2*m[i][j]);
    t = sign(theta) / (fabs(theta) + sqrt(pow(theta, 2) + 1));
    return t;
}

double sign(double x) {
    /* given double returns it's sign. Note: sign(0) = 1.*/
    double s;
    if(x == 0){
        return 1;
    }
    s = (double)(x > 0) - (x < 0);
    return s;
}

int *max_off_diagonal(double **m){
/**
 * Finds the indices of the maximum double in a matrix
 *
 * @param m - n*n matrix
 * @return [i,j] - 2-array of int representing row and column of the maximal double
 */
    int *res; double max = 0; int i,j;
    res =  (int*) calloc(2, sizeof(int));
    /* Just sanity check if the matrix is already all zeroes:*/
    res[0] = 0;
    res[1] = 1;
    if (res == NULL) {
        print_error();
    }
    for(i=0;i<n;i++){
        for(j=i+1;j<n;j++){
            if (fabs(m[i][j]) > fabs(max)){
                max = m[i][j];
                res[0] = i;
                res[1] = j;
            }
        }
    }
    return res;
}

void calc_a_tag(double **a, double **a_tag, double s, double c, int i, int j){
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
    a_tag[j][i] = 0;
}

void calc_p(double **p, int i, int j, double s, double c){
/**
 * Calculates Rotation Matrices in place using heuristics
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
        p[r][j] = p[r][j] + s*(temp - tau*p[r][j]);
    }
}

void copy_matrix(double **m, double **c) {
    /* Deep copy a n*n matrix*/
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            c[i][j] = m[i][j];
        }
    }
}

double **identity_matrix(int k) {
    /* return an N*N identity matrix*/
    double **identity; int i;
    identity = allocate_data(n, n);
    for (i = 0; i < k; i++) {
        identity[i][i] = 1;
    }
    return identity;
}

double *get_diagonal(double **m) {
    /* return an N*N identity matrix*/
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

double **degree_to_diagonal_matrix(double *degree){
    /* Use to convert the 1D array to 2D array*/
    double **diagonal_degree_matrix_res; int i;
    diagonal_degree_matrix_res = allocate_data(n, n);
    for (i=0 ; i<n; i++){
            diagonal_degree_matrix_res[i][i] = degree[i];
    }
    return diagonal_degree_matrix_res;
}


int idx_cmp( const void* a, const void* b)
{
    /* Sort indexed double array*/
    indexed_double *x = (indexed_double *) a;
    indexed_double *y = (indexed_double *) b;

    if ( (*x).v == (*y).v ) {
        return 0;
    }
    else if ( (*x).v < (*y).v )
    {return -1;}
    return 1;
}


int eigengap_hueuristic(const double *eigenvaleus){
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
void sort_eigenvalues_and_vectors(const double *eigenvalues, double **eigenvectors,
                                         double * s_eigenvalues, double ** s_eigenvectors){
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
    for(i = 0; i<n; i++){
        s_eigenvalues[i] = eigenvalues[eigenvalues_idx[i].i];
        s_eigenvectors[i] = eigenvectors[eigenvalues_idx[i].i];
    }
    free(eigenvalues_idx);
}

double **calculate_T(double **eigenvectors, int k){
    double **T, *norms, sum;
    int i,j;
    T = allocate_data(n, k);
    norms = (double *) calloc(k, sizeof (double));
    if (norms == NULL) {
        print_error();
    }
    /* Calculate norms:*/
    for (j=0;j<k;j++){
        sum = 0;
        for(i=0;i<n;i++){
            sum+= pow(eigenvectors[i][j], 2);
        }
        norms[j] = sqrt(sum);
    }
    for (i=0;i<n;i++){
        for(j=0;j<k;j++){
            T[i][j] = (eigenvectors[i][j] / norms[j]);
        }
    }
    free(norms);
    return T;
}


double **spectral_clustrering(double **datapoints){
/**
 * Preform the full spectral k-means algorithm, return optimal k.
 *
 * @param datapoints - 2D array of datapoints size n*d.
 * @invariant - n and d are initialized.
  *@result k - optimal k for k-means algorithm.
 */    
    double **weights, *degree, **laplacian, *eigenvalues, **eigenvectors,
            *s_eigenvalues, **s_eigenvectors, **T;
    int k;
    s_eigenvalues = (double *) calloc(n, sizeof (double));
    if (s_eigenvalues == NULL) {
        print_error();
    }
    s_eigenvectors = allocate_data(n, n);
    weights = weight_adj_matrix(datapoints);
    degree = diagonal_degree_matrix(weights, 1);
    laplacian = normalized_laplacian(weights, degree);
    eigenvectors = jacobi_function(laplacian, EPSILON);
    eigenvalues = get_diagonal(laplacian);
    sort_eigenvalues_and_vectors(eigenvalues, eigenvectors, s_eigenvalues, s_eigenvectors);
    k = eigengap_hueuristic(s_eigenvalues);
    T = calculate_T(s_eigenvectors, k);
    free(s_eigenvalues);
    free(s_eigenvectors);
    free(weights);
    free(degree);
    free(laplacian);
    free(eigenvectors);
    free(eigenvalues);
    return T;
}


/*kmeans from first and second exc:*/
double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon) {
    double max_change; int iterations=0;
    do{
        max_change = update_centroids(centroids, datapoints, k);
        iterations++;
    } while ((max_change >= epsilon) && (iterations < max_iter));
    return centroids;
}

double update_centroids(double **centroids, double **datapoints, int k){
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
        update_cumulative_sums(datapoints[i], cumulative_sums[chosen_m_idx]);
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

void update_cumulative_sums(double *arr, double *cumulative_sum){
    int i;
    for (i=0; i<d; i++){
        cumulative_sum[i] += arr[i];
    }
}


int main(int argc, char *argv[]) {
    char *file_name, *goal;
    double **eigen_vectors, **diagonal_degree_matrix_res, **lap_res, *degree,
        **datapoints, **weight_adj_matrix_res, *eigenvalues;
    
    if (argc == 3) {
        goal = argv[1];
        file_name = argv[2];
    }
    else{
        print_invalid_input();
    }
    datapoints = parse_file(file_name);
    if(strcmp(goal, "wam") == 0) {
        weight_adj_matrix_res = weight_adj_matrix(datapoints);
        print_matrix(weight_adj_matrix_res, n, n);
    }
    else if(strcmp(goal, "ddg") == 0) {

        weight_adj_matrix_res = weight_adj_matrix(datapoints);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 0);
        diagonal_degree_matrix_res = degree_to_diagonal_matrix(degree);
        print_matrix(diagonal_degree_matrix_res, n, n);
    }
    else if(strcmp(goal, "lnorm") == 0) {
        weight_adj_matrix_res = weight_adj_matrix(datapoints);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 1);
        lap_res = normalized_laplacian(weight_adj_matrix_res, degree);
        print_matrix(lap_res, n, n);
    }
    else if(strcmp(goal, "jacobi") == 0) {
        eigen_vectors = jacobi_function(datapoints, EPSILON);
        eigenvalues = get_diagonal(datapoints);
        print_arr(eigenvalues, n);
        print_matrix(eigen_vectors, n, n);
    }
   /* TODO: FOR tests only - delete before submission*/
    else if(strcmp(goal, "spk") == 0) {
        eigen_vectors = spectral_clustrering(datapoints);
        print_matrix(eigen_vectors, n, n);
    }
    else{
            printf("Invalid Input!");
            return 1;
    }
    return 0;
}

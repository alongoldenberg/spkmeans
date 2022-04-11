#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "utils.c"
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))


int d; int n;

void print_error();
double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon);
double update_centroids(double **centroids, double **datapoints, int k, int n, int d);
void update_cumulative_sums(double *arr, double *cumulative_sum);

void print_error(){
    printf("An Error Has Occurred!\n");
    exit(1);
}

double **allocate_data(int n, int d) {
    double *p; double **a; int i;
    p = (double*)calloc(n*d, sizeof(double));
    a = (double**) calloc(n,sizeof(double *));
    if(a == NULL || p == NULL){
        print_error();
    }
    for(i=0; i<n; i++) {
        a[i] = p + i*d;
    }
    return a;
}


double **kmeans(double **datapoints, double **centroids, int k, int max_iter, double epsilon) {
    double max_change; int iterations=0;
    do{
        max_change = update_centroids(centroids, datapoints, k, n , d);
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
            if (distance(centroids[j], datapoints[i], d) <
                distance(centroids[chosen_m_idx], datapoints[i], d)){
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

/*
 * Below - Python C-API implementation:
 */

void initalize_datapoints(PyObject *dp_pointer, double **datapoints){
    PyObject *t;
    int i, j;
    for(i=0; i<n;i++){
        t = PyList_GetItem(dp_pointer, i);
        for(j=0;j<d;j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(t, j));

        }
    }
}

void initalize_centroids(PyObject *centroids_pointer, double **centroids, int k){
    PyObject *t; int i, j;
    for(i=0; i<k;i++){
        t = PyList_GetItem(centroids_pointer, i);
        for(j=0;j<d;j++)
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(t, j));
    }
}

void put_centroids(PyObject *py_centroids, double **centroids, int k){
    PyObject *t; int i, j;
    for(i=0; i<k;i++){
        t = PyList_New(d);
        if(t==NULL){print_error();}
        for(j=0;j<d;j++) {
            PyList_SetItem(t, j, PyFloat_FromDouble(centroids[i][j]));
        }
        PyList_SetItem(py_centroids, i, t);
    }
}

PyObject* mykmeanssp(PyObject *self, PyObject *args){
    PyObject *dp_pointer, *centroids_pointer, *py_centroids;
    double **datapoints, **centroids, epsilon;
    int k, max_iter;
    /* Validate input: */
    if(!PyArg_ParseTuple(args, "OOiidii", &dp_pointer, &centroids_pointer, &k, &max_iter, &epsilon, &n, &d)) {
        return NULL;
    }

    /* Initalize arrays and run kmeans */
    datapoints = allocate_data(n, d);
    centroids = allocate_data(k, d);
    initalize_datapoints(dp_pointer, datapoints);
    initalize_centroids(centroids_pointer, centroids, k);
    centroids = kmeans(datapoints, centroids, k, max_iter, epsilon);

    /* handle returning result as list */

    py_centroids = PyList_New(k);
    if (py_centroids == NULL)
        print_error();
    put_centroids(py_centroids, centroids, k);
    free(centroids[0]);
    free(centroids);
    free(datapoints[0]);
    free(datapoints);
    return py_centroids;
}

PyMethodDef kmeans_methods[] = {
        {"fit",
                (PyCFunction) mykmeanssp,
                     METH_VARARGS,
                            PyDoc_STR("Runs kmeans")},
        {NULL, NULL, 0, NULL}
};

struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        kmeans_methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    return PyModule_Create(&_moduledef);
}

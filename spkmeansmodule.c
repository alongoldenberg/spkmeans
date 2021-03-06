#define PY_SSIZE_T_CLEAN  
#include <Python.h>                                  
#include <math.h>         
#include "utils.h"
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static double** data_py_to_c_type(PyObject* datapoints_py_type,
                                 int data_size_n, int data_dimension);
static PyObject* data_c_to_py_type(double** data_c_type,
                                   int data_size_n, int data_dimension);
static PyObject* get_goal_capi(PyObject *self, PyObject *args);
static PyObject* calc_kmeans_capi(PyObject *self, PyObject *args);
static int calculate_k(double **datapoints, int n, int d);


static int calculate_k(double **datapoints, int n, int d){
/**
 * Find k using eigengap_hueuristic function.
 * 
 * @param datapoints - 2D array of datapoints size n*d.
 * @invariant - n and d are initialized.
  *@result k
 */
    double **weights, *degree, **laplacian, *eigenvalues, **eigenvectors,
            *s_eigenvalues, **s_eigenvectors;
    int k;
    s_eigenvalues = (double *) calloc(n, sizeof (double));
    if (s_eigenvalues == NULL) {
        print_error();
    }
    s_eigenvectors = allocate_data(n, n);
    weights = weight_adj_matrix(datapoints, n, d);
    degree = diagonal_degree_matrix(weights, 1, n);
    laplacian = normalized_laplacian(weights, degree, n);
    eigenvectors = jacobi_function(laplacian, EPSILON, n);
    eigenvalues = get_diagonal(laplacian, n);
    sort_eigenvalues_and_vectors(eigenvalues, eigenvectors,
                                s_eigenvalues, s_eigenvectors, n);
    k = eigengap_hueuristic(s_eigenvalues, n);
    
    free(weights[0]);
    free(weights);
    free(laplacian[0]);
    free(laplacian);
    free(degree);
    free(eigenvectors[0]);
    free(eigenvectors);
    free(s_eigenvectors);
    free(eigenvalues);
    free(s_eigenvalues);
    return k;
}


static double** data_py_to_c_type(PyObject* datapoints_py_type,
                                  int data_size_n, int data_dimension){
/**
 * Create C type matrix and copy data to it from input
 * @param datapoints_py_type - 2D array of datapoints size data_size_n*data_dimension.
 * @invariant - data_size_n and data_dimension are initialized.
 * @result c_data_pointer - requested matrix
 * 
 */
    double **c_data_pointer; Py_ssize_t i, j;
    PyObject *py_point_vector;
    c_data_pointer = allocate_data(data_size_n, data_dimension);
    for (i = 0; i < data_size_n; i++) {
        py_point_vector = PyList_GetItem(datapoints_py_type, i);
        if (!PyList_Check(py_point_vector)) {
            free(c_data_pointer[0]);  
            free(c_data_pointer);
            return NULL;
            }
        for (j=0; j < data_dimension; j++){
            c_data_pointer[i][j] = PyFloat_AsDouble(PyList_GetItem(
                                                    py_point_vector, j));
        }
    }
    return c_data_pointer;
}


static PyObject* data_c_to_py_type(double** data_c_type,
                                   int data_size_n, int data_dimension){
/**
 * Create Py type matrix and copy data to it from input
 * @param data_c_type - 2D array of datapoints size data_size_n*data_dimension.
 * @invariant - data_size_n and data_dimension are initialized.
 * @result py_data_pointer - requested matrix
 * 
 */
    PyObject* py_data_pointer = PyList_New(0);
    int i, j;
    PyObject* py_point_vector;
    for (i = 0; i < data_size_n; i++) {
        py_point_vector = PyList_New(0);
        for (j=0; j< data_dimension; j++){
            PyList_Append(py_point_vector, PyFloat_FromDouble(data_c_type[i][j]));
        }
        PyList_Append(py_data_pointer, py_point_vector);
    }
    return py_data_pointer;
}


static PyObject* get_goal_capi(PyObject *self, PyObject *args){
/**
 * Get goal set by the user through spkmeans.py (handle each goal separately).
 * return PyObject of the results
 */
    PyObject *datapoints_py_type, *result;
    char* my_goal_c_type;
    int n, d, k; 
    double **datapoints; 
    /* parse arguments: */
    if(!PyArg_ParseTuple(args, "iiisO:get_goal", &n, &d, &k,
                         &my_goal_c_type, &datapoints_py_type)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        return NULL;
    }
    /* create c type object of datapoints: */
    datapoints = data_py_to_c_type(datapoints_py_type, n, d);   
    
    /* handle goals: */
    if  (strcmp("spk", my_goal_c_type)==0){
        double **T;
        T = spectral_clustrering(datapoints, n, d, k);
        if(k==0){
            k = calculate_k(datapoints, n,d);
        }
        if(k==1){
            print_error();
        }
        result = data_c_to_py_type(T, n, k);
        free(datapoints[0]);
        free(datapoints);
        free(T[0]);
        free(T);
        return Py_BuildValue("O", result);
    }
    else if(strcmp(my_goal_c_type, "wam") == 0) {
        double **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        result = data_c_to_py_type(weight_adj_matrix_res, n, n);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
    }
    else if(strcmp(my_goal_c_type, "ddg") == 0) {
        double **diagonal_degree_matrix_res, *degree, **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 0, n);
        diagonal_degree_matrix_res = degree_to_diagonal_matrix(degree, n);
        result = data_c_to_py_type(diagonal_degree_matrix_res, n, n);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
        free(degree);
        free(diagonal_degree_matrix_res[0]);
        free(diagonal_degree_matrix_res);
    }
    else if(strcmp(my_goal_c_type, "lnorm") == 0) {
        double **lap_res, *degree, **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 1, n);
        lap_res = normalized_laplacian(weight_adj_matrix_res, degree, n);
        result = data_c_to_py_type(lap_res, n, n);
        free(degree);
        free(weight_adj_matrix_res[0]);
        free(weight_adj_matrix_res);
        free(lap_res[0]);
        free(lap_res);
    }
    else if(strcmp(my_goal_c_type, "jacobi") == 0) {
        double **eigen_vectors;
        PyObject *eigen_vectors_py, *datapoints_py_type; 

        /*here datapoints is symitric matrix with eignvalues over the diagonal*/
        eigen_vectors = jacobi_function(datapoints, EPSILON, n); 
        eigen_vectors_py = data_c_to_py_type(eigen_vectors, n, n);
        datapoints_py_type = data_c_to_py_type(datapoints, n, n);
        free(eigen_vectors[0]);
        free(eigen_vectors);
        free(datapoints[0]);
        free(datapoints);
        /* in this case only the function returns two matrixes */
        return Py_BuildValue("OO", datapoints_py_type, eigen_vectors_py);
    }
    else{ /* if the goal is not valid */
        printf("Invalid Input!");
        free(datapoints[0]);
        free(datapoints);
        return NULL;
    }
    free(datapoints[0]);
    free(datapoints);
    return Py_BuildValue("O", result);
}


static PyObject* calc_kmeans_capi(PyObject *self, PyObject *args){
/**
 * update centroids according to kmeans algorithm.
 * this function parse PyObjects to c objectes, calculate centroids,
 * and return the new centroids as PyObjects
 */
    PyObject* datapoints_py_type;
    PyObject* centroids_py_type;
    double** datapoints;
    double** centroids;
    int n, d, k, eps;
    eps = 0; // Default Value
    if(!PyArg_ParseTuple(args, "iiiOO:calc_kmeans", &n, &d, &k,
                         &datapoints_py_type, &centroids_py_type)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        return NULL;
    }
    if (!PyList_Check(centroids_py_type)){
        return NULL;
    }
    datapoints = data_py_to_c_type(datapoints_py_type, n, d);
    centroids = data_py_to_c_type(centroids_py_type, k, d);
    centroids = kmeans(datapoints, centroids, k, MAX_ITER, eps, n, d);
    centroids_py_type = data_c_to_py_type(centroids, k, d);
    free(centroids[0]);
    free(centroids);
    free(datapoints[0]);
    free(datapoints);
    return Py_BuildValue("O", centroids_py_type);
}


/*
 * The methods this module has.
 * We will use it in the next structure
 */
static PyMethodDef spkmethods[] = {
        {"get_goal",
                (PyCFunction) get_goal_capi,
                     METH_VARARGS,
                        PyDoc_STR("calculate goal set by py program")},
        {"kmeans",
                (PyCFunction) calc_kmeans_capi,
                     METH_VARARGS,
                        PyDoc_STR("return centroids using kmeans algorithm")},
        {NULL, NULL, 0, NULL}
};

/* Initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans", /* name of module */
        NULL,
        -1,
        spkmethods
};


PyMODINIT_FUNC
PyInit_myspkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}



#define PY_SSIZE_T_CLEAN
#define PY_SSIZE_T_CLEAN  /* For all # variants of unit formats (s#, y#, etc.) use Py_ssize_t rather than int. */
#include <Python.h>       /* MUST include <Python.h>, this implies inclusion of the following standard headers:
                             <stdio.h>, <string.h>, <errno.h>, <limits.h>, <assert.h> and <stdlib.h> (if available). */
#include <math.h>         /* include <Python.h> has to be before any standard headers are included */

#include "utils.h"
#include "spkmeans.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

static int calculate_k(double **datapoints, int n, int d);
static double** data_py_to_c_type(PyObject* datapoints_py_type, int data_size_n, int data_dimension);
static PyObject* data_c_to_py_type(double** data_c_type, int data_size_n, int data_dimension);
static PyObject* get_goal_capi(PyObject *self, PyObject *args);
static PyObject* calc_kmeans_capi(PyObject *self, PyObject *args);




static int calculate_k(double **datapoints, int n, int d){
/**
 * Preform the full spectral k-means algorithm, return k first eigenvectors.
 *
 * @param datapoints - 2D array of datapoints size n*d.
 * @invariant - n and d are initialized.
  *@result T - requested matrix 
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
    sort_eigenvalues_and_vectors(eigenvalues, eigenvectors, s_eigenvalues, s_eigenvectors, n);
    k = eigengap_hueuristic(s_eigenvalues, n);
    return k;
}



static double** data_py_to_c_type(PyObject* datapoints_py_type, int data_size_n, int data_dimension){
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
            c_data_pointer[i][j] = PyFloat_AsDouble(PyList_GetItem(py_point_vector, j));
        }
    }
    return c_data_pointer;
}

static PyObject* data_c_to_py_type(double** data_c_type, int data_size_n, int data_dimension){
    PyObject* py_data_pointer = PyList_New(0);
    int i, j;
    PyObject* py_point_vector;
    for (i = 0; i < data_size_n; i++) {
        py_point_vector = PyList_New(0);
        for (j=0; j< data_dimension; j++){
            PyList_Append(py_point_vector, PyFloat_FromDouble(data_c_type[j][i]));
        }
        PyList_Append(py_data_pointer, py_point_vector);
    }
    return py_data_pointer;
}


static PyObject* get_goal_capi(PyObject *self, PyObject *args){
    PyObject *datapoints_py_type, *result;
    char* my_goal_c_type;
    int n, d;
    double **datapoints; 
    if(!PyArg_ParseTuple(args, "iisO:get_goal", &n, &d, &my_goal_c_type, &datapoints_py_type)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        return NULL;
    }
       
    datapoints = data_py_to_c_type(datapoints_py_type, n, d);   
    
    if  (strcmp("spk_T_and_k", my_goal_c_type)==0){
        double **T;
        int k;
        T = spectral_clustrering(datapoints, n, d);
        k = calculate_k(datapoints,  n, d);
        // print_matrix(eigen_vectors, n, d);
        result = data_c_to_py_type(T, n, k);
        return Py_BuildValue("Oi", result, k);
    }
    else if(strcmp(my_goal_c_type, "wam") == 0) {
        double **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        //print_matrix(weight_adj_matrix_res, n, n);
        result = data_c_to_py_type(weight_adj_matrix_res, n, n);
    }
    else if(strcmp(my_goal_c_type, "ddg") == 0) {
        double **diagonal_degree_matrix_res, *degree, **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 0, n);
        diagonal_degree_matrix_res = degree_to_diagonal_matrix(degree, n);
        //print_matrix(diagonal_degree_matrix_res, n, n);
        result = data_c_to_py_type(diagonal_degree_matrix_res, n, n);
    }
    else if(strcmp(my_goal_c_type, "lnorm") == 0) {
        double **lap_res, *degree, **weight_adj_matrix_res;
        weight_adj_matrix_res = weight_adj_matrix(datapoints, n, d);
        degree = diagonal_degree_matrix(weight_adj_matrix_res, 1, n);
        lap_res = normalized_laplacian(weight_adj_matrix_res, degree, n);
        //print_matrix(lap_res, n, n);
        result = data_c_to_py_type(lap_res, n, n);
    }
    else if(strcmp(my_goal_c_type, "jacobi") == 0) {
        double **eigen_vectors;
        PyObject *eigen_vectors_py;
        eigen_vectors = jacobi_function(datapoints, EPSILON, n);  /*here datapoints is symitric matrix with eignvalues over the diagonal*/
        // print_matrix(datapoints, n, n);
        eigen_vectors_py = data_c_to_py_type(eigen_vectors, n, n);
        return Py_BuildValue("OO", data_c_to_py_type(datapoints, n, n), eigen_vectors_py);
    }
    else{
        printf("Invalid Input!");
        return NULL;
    }
    return Py_BuildValue("O", result);
}


static PyObject* calc_kmeans_capi(PyObject *self, PyObject *args){
    PyObject* datapoints_py_type;
    PyObject* centroids_py_type;
    double** datapoints;
    double** centroids;
    int n, d, k;
    if(!PyArg_ParseTuple(args, "iiiOO:calc_kmeans", &n, &d, &k, &datapoints_py_type, &centroids_py_type)) {
        puts("problem 1");
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        puts("problem 2");
        return NULL;
    }
    if (!PyList_Check(centroids_py_type)){
        puts("problem 3");
        return NULL;
    }
        
    datapoints = data_py_to_c_type(datapoints_py_type, n, d);
    centroids = data_py_to_c_type(centroids_py_type, k, d);
    centroids = kmeans(datapoints, centroids, k, MAX_ITER, EPSILON, n, d);
    centroids_py_type = data_c_to_py_type(centroids, k, d);
    return Py_BuildValue("O", centroids_py_type);
}


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

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
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

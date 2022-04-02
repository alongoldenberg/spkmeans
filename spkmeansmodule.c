//
// Created by Alon Goldenberg & Ofir Nissan.
//

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "utils.c"
#include <spkmeans.h>


static double** data_py_to_c_type(PyObject* datapoints_py_type, int data_size_n, int data_dimension){
    **double c_data_pointer; int i, j;
    PyObject* py_point_vector;
    c_data_pointer = allocate_data(data_size_n, data_dimension);
    for (i = 0; i < data_size_n; i++) {
        py_point_vector = PyList_GetItem(datapoints_py_type, i);
        if (!PyList_Check(row)) {
            free(c_data_pointer[0]);  // free up the memory before leaving
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
        py_point_vector = PyList_New(0)
        for (j=0; j< data_dimension; j++){
            PyList_Append(py_point_vector, PyFloat_FromDouble(data_c_type[i][j]));
        }
        PyList_Append(py_data_pointer, py_point_vector);
        printf("check first coordinate of first vector in py list %d", py_data_pointer[0][0]); // sanity check
    }
    return py_data_pointer;
}


static PyObject* get_goal_capi(PyObject *self, PyObject *args){
    PyObject* datapoints_py_type;
    double** datapoints;
    char* my_goal_c_type;
    int data_size_n, data_dimension, k;
    double** result, weight_adj_matrix_res, lap, T;
    double* degree;
    if(!PyArg_ParseTuple(args, "iisO:get_goal", &data_size_n, &data_dimension, &my_goal_c_type, &datapoints_py_type)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        return NULL;
    }
    //initiate global variables d and n at spkmeans.c:
    init_d_and_n(data_size_n, data_dimension);
    //convert datapoints:
    datapoints = data_py_to_c_type(datapoints_py_type, data_size_n, data_dimension);
    //calculations using spkmeans.c
    weight_adj_matrix_res = weight_adj_matrix(datapoints);
    degree = diagonal_degree_matrix(weight_adj_matrix_res);
    lap = normalized_laplacian(weight_adj_matrix_res, degree);
    T = find_T(datapoints);
    k = calculate_k(datapoints);
    //cases:
    switch (my_goal_c_type)
    {
        case "spk_T_and_k":
            result = T;
            return Py_BuildValue("Oi", result, k);
        case "wam":
            result = weight_adj_matrix_res;
            break;
        case "ddg":
            result = degree_to_diagonal_matrix(degree);
            break;
        case "lnorm":
            result = lap;
            break;
        case "jacobi":
            result = jacobi_function(datapoints, EPSILON);
        default:
            printf("Invalid Input!");
            return NULL
    }
    return Py_BuildValue("O", result);
}


static PyObject* calc_kmeans_capi(PyObject *self, PyObject *args){
    PyObject* datapoints_py_type;
    PyObject* centroids_py_type;
    double** datapoints;
    double** centroids;
    int data_size_n, data_dimension, k;
    double** result, weight_adj_matrix_res, lap;
    if(!PyArg_ParseTuple(args, "iiiOO:get_goal", &data_size_n, &data_dimension, &k, &datapoints_py_type, &centroids_py_type)) {
        return NULL; 
    }
    if (!PyList_Check(datapoints_py_type)){
        return NULL;
    }
    if (!PyList_Check(centroids_py_type)){
        return NULL;
    }
    //initiate global variables d and n at spkmeans.c:
    init_d_and_n(data_size_n, data_dimension);
    //create c_type data and centroids from python type:
    data = data_py_to_c_type(datapoints_py_type, data_size_n, data_dimension);
    centroids = data_py_to_c_type(centroids_py_type, k, data_dimension);
    //calc centroids and return py object:
    centroids = kmeans(datapoints, centroids, k, MAX_ITTER, EPSILON);
    centroids_py_type = data_c_to_py_type(centroids, k, data_dimension);
    return Py_BuildValue("O", centroids_py_type);
}


static PyMethodDef spk_methods[] = {
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

static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        spk_methods
};

PyMODINIT_FUNC PyInit_myspkmeans(void)
{
    return PyModule_Create(&_moduledef);
}
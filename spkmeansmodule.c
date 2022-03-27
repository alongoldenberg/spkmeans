//
// Created by Alon Goldenberg & Ofir Nissan.
//

#include <Python.h>
#include "utils.c"

static void initalize_datapoints(PyObject *dp_pointer, double **datapoints){
    PyObject *t;
    int i, j;
    for(i=0; i<n;i++){
        t = PyList_GetItem(dp_pointer, i);
        for(j=0;j<d;j++) {
            datapoints[i][j] = PyFloat_AsDouble(PyList_GetItem(t, j));

        }
    }
}


static PyObject* myspkmeans(PyObject *self, PyObject *args){
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
    //put_centroids(py_centroids, centroids, k);
    free(centroids[0]);
    free(centroids);
    free(datapoints[0]);
    free(datapoints);
    return py_centroids;
}

static PyMethodDef kmeans_methods[] = {
        {"fit",
                (PyCFunction) myspkmeans,
                     METH_VARARGS,
                        PyDoc_STR("Runs kmeans")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef _moduledef = {
        PyModuleDef_HEAD_INIT,
        "myspkmeans",
        NULL,
        -1,
        kmeans_methods
};

PyMODINIT_FUNC PyInit_myspkmeans(void)
{
    return PyModule_Create(&_moduledef);
}
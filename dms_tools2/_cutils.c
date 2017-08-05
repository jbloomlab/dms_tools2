// Fast C versions of some functions in dms_tools2.utils
// Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static PyObject *
buildReadConsensus(PyObject *self, PyObject *args)
{
    const char *command;

    if (! PyArg_ParseTuple(args, "s", &command))
        return NULL;

    PyErr_SetString(PyExc_RuntimeError, "not yet implemented");
    return NULL;
}


static PyMethodDef cutilsMethods[] = {
    {"buildReadConsensus", buildReadConsensus, METH_VARARGS,
            "Fast version of `dms_tools2.utils.buildReadConsensus`."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef cutilsmodule = {
    PyModuleDef_HEAD_INIT,
    "_cutils",
    "Fast implementations of some functions in `dms_tools2.utils`.",
    -1,
    cutilsMethods
};

PyMODINIT_FUNC
PyInit__cutils(void)
{
    return PyModule_Create(&cutilsmodule);
}

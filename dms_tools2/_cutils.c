// Fast C versions of some functions in dms_tools2.utils
// Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static PyObject *
lowQtoN(PyObject *self, PyObject *args)
{
    // define variables
    PyObject *py_newr;
    int minq;
    const char *r, *q;
    Py_ssize_t rlen, i;

    // parse arguments
    if (! PyArg_ParseTuple(args, "ssC", &r, &q, &minq)) {
        return NULL;
    }
    rlen = strlen(r);
    if (rlen != strlen(q)) {
        PyErr_SetString(PyExc_ValueError, "r and q not of same length");
        return NULL;
    }

    // build up new string
    char *newr = PyMem_New(char, rlen + 1);
    if (newr == NULL) {
        PyErr_SetString(PyExc_MemoryError, "cannot allocate newr");
        return NULL;
    }
    newr[rlen] = '\0'; // string termination character
    for (i = 0; i < rlen; i++) {
        if (q[i] >= minq) {
            newr[i] = r[i];
        } else {
            newr[i] = 'N';
        }
    }
    py_newr = PyUnicode_FromString(newr);    
    PyMem_Del(newr);
    return py_newr;
}

static PyObject *
buildReadConsensus(PyObject *self, PyObject *args)
{
    // define variables
    PyObject *r, *py_consensus;
    PyListObject *reads;
    long minreads;
    double minconcur, mincount;
    Py_ssize_t nreads, iread, maxrlen, rlen, i;
    const char *rchar;

    // arrays to hold counts
    const int maxallowrlen = 2000; // allow reads up to this length
    long countsA[maxallowrlen], countsC[maxallowrlen];
    long countsG[maxallowrlen], countsT[maxallowrlen];
    long countstot[maxallowrlen];

    // parse arguments
    if (! PyArg_ParseTuple(args, "O!ld", &PyList_Type, &reads, 
            &minreads, &minconcur)) {
        return NULL;
    }
    nreads = PyList_GET_SIZE(reads);
    if (nreads < 1) {
        PyErr_SetString(PyExc_ValueError, "reads has no reads");
        return NULL;
    }
    maxrlen = 0;
    for (iread = 0; iread < nreads; iread++) {
        r = PyList_GET_ITEM(reads, iread);
        if (! PyUnicode_Check(r)) {
            PyErr_SetString(PyExc_ValueError, "entry in reads not unicode");
            return NULL;
        }
        rlen = PyUnicode_GET_LENGTH(r);
        if (rlen > maxrlen) {
            maxrlen = rlen;
        }
    }

    // Count nucleotide occurrences
    for (i = 0; i < maxrlen; i++) {
        countsA[i] = 0;
        countsC[i] = 0;
        countsG[i] = 0;
        countsT[i] = 0;
        countstot[i] = 0;
    }
    for (iread = 0; iread < nreads; iread++) {
        rchar = PyUnicode_AsUTF8AndSize(PyList_GET_ITEM(reads, iread), &rlen);
        for (i = 0; i < rlen; i++) {
            countstot[i]++;
            switch (rchar[i]) {
                case 'A' : countsA[i]++;
                           break;
                case 'C' : countsC[i]++;
                           break;
                case 'G' : countsG[i]++;
                           break;
                case 'T' : countsT[i]++;
                           break;
                case 'N' : countstot[i]--;
                           break;
                default : PyErr_SetString(PyExc_ValueError, "invalid nt");
                          return NULL;
            }
        }
    }
        
    // build consensus sequence
    char *rconsensus = PyMem_New(char, maxrlen + 1);
    if (rconsensus == NULL) {
        PyErr_SetString(PyExc_MemoryError, "cannot allocate rconsensus");
        return NULL;
    }
    rconsensus[maxrlen] = '\0'; // string termination character
    for (i = 0; i < maxrlen; i++) {
        mincount = minconcur * countstot[i];
        if (mincount < minreads) {
            mincount = minreads;
        }
        if (countsA[i] >= mincount) {
            rconsensus[i] = 'A';
        } else if (countsC[i] >= mincount) {
            rconsensus[i] = 'C';
        } else if (countsG[i] >= mincount) {
            rconsensus[i] = 'G';
        } else if (countsT[i] >= mincount) {
            rconsensus[i] = 'T';
        } else {
            rconsensus[i] = 'N';
        }
    }
    py_consensus = PyUnicode_FromString(rconsensus);    
    PyMem_Del(rconsensus);
    return py_consensus;
}


static PyMethodDef cutilsMethods[] = {
    {"buildReadConsensus", buildReadConsensus, METH_VARARGS,
            "Fast version of `dms_tools2.utils.buildReadConsensus`."},
    {"lowQtoN", lowQtoN, METH_VARARGS,
            "Fast version of `dms_tools2.utils.lowQtoN`."},
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

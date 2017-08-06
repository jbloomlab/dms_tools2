// Fast C versions of some functions in dms_tools2.utils
// Written by Jesse Bloom.
//
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


static PyObject *
alignSubamplicon(PyObject *self, PyObject *args)
{
    // define variables
    const char *refseq, *r1, *r2, *chartype;
    double maxmuts, maxN;
    long refseqstart, refseqend, i, j;
    long startcodon, nmuts, codonshift;
    char mutnt;
    int hasN, hasmut;
    PyObject *py_subamplicon;

    // parse arguments
    if (! PyArg_ParseTuple(args, "ssslldds", &refseq, &r1, &r2,
            &refseqstart, &refseqend, &maxmuts, &maxN, &chartype)) {
        return NULL;
    }
    if (strcmp(chartype, "codon")) {
        PyErr_SetString(PyExc_ValueError, "chartype not codon");
        return NULL;
    }
    long len_subamplicon = refseqend - refseqstart + 1;
    long len_r1 = strlen(r1);
    long len_r2 = strlen(r2);
    long len_subamplicon_minus_len_r2 = len_subamplicon - len_r2;

    // build subamplicon
    char *subamplicon = PyMem_New(char, len_subamplicon + 1);
    if (subamplicon == NULL) {
        PyErr_SetString(PyExc_MemoryError, "cannot allocate subamplicon");
        return NULL;
    }
    subamplicon[len_subamplicon] = '\0'; // string termination character
    long nN = 0;
    for (i = 0; i < len_subamplicon; i++) {
        if (i < len_subamplicon_minus_len_r2) { // site not in r2
            if (i < len_r1) { // site in r1
                subamplicon[i] = r1[i];
            } else { // site not in r1
                subamplicon[i] = 'N';
            }
        } else { // site in r2
            if (i < len_r1) { // site in r1
                if (r1[i] == r2[i - len_subamplicon_minus_len_r2]) {
                    subamplicon[i] = r1[i];
                } else if (r1[i] == 'N') {
                    subamplicon[i] = r2[i - len_subamplicon_minus_len_r2];
                } else if (r2[i - len_subamplicon_minus_len_r2] == 'N') {
                    subamplicon[i] = r1[i];
                } else {
                    subamplicon[i] = 'N';
                }
            } else { // site not in r1
               subamplicon[i] = r2[i - len_subamplicon_minus_len_r2]; 
            }
        }
        if (subamplicon[i] == 'N') {
            nN++;
            if (nN > maxN) {
                Py_RETURN_FALSE;
            }
        }
    }

    // look for excessive mutations
    if (! strcmp(chartype, "codon")) {
        switch (refseqstart % 3) {
            case 1 : startcodon = (refseqstart + 2) / 3;
                     codonshift = 0;
                     break;
            case 2 : startcodon = (refseqstart + 1) / 3 + 1;
                     codonshift = 2;
                     break;
            case 0 : startcodon = refseqstart / 3 + 1;
                     codonshift = 1;
                     break;
            default : PyErr_SetString(PyExc_ValueError, "invalid case");
                      return NULL;
        }
        nmuts = 0;
        for (i = startcodon; i < (refseqend / 3 + 1); i++) {
            hasN = 0;
            hasmut = 0;
            for (j = 0; j < 3; j++) {
                mutnt = subamplicon[3 * (i - startcodon) + codonshift + j];
                if (mutnt == 'N') {
                    hasN = 1;
                    break;
                } else if (mutnt != refseq[3 * i - 3 + j]) {
                    hasmut = 1;
                }
            }
            if (hasmut && (! hasN)) {
                nmuts++;
                if (nmuts > maxmuts) {
                    Py_RETURN_FALSE;
                }
            }
        }
    } else {
        PyErr_SetString(PyExc_ValueError, "invalid chartype");
        return NULL;
    }

    // return subamplicon
    py_subamplicon = PyUnicode_FromString(subamplicon);
    PyMem_Del(subamplicon);
    return py_subamplicon;
}


static PyObject *
reverseComplement(PyObject *self, PyObject *args)
{
    // define variables
    PyObject *py_rc;
    const char *s;
    size_t slen, i;

    // parse arguments
    if (! PyArg_ParseTuple(args, "s", &s)) {
        return NULL;
    }
    slen = strlen(s);

    // build up new string
    char *rc = PyMem_New(char, slen + 1);
    if (rc == NULL) {
        PyErr_SetString(PyExc_MemoryError, "cannot allocate rc");
        return NULL;
    }
    rc[slen] = '\0'; // string termination character
    for (i = 0; i < slen; i++) {
        switch (s[slen - 1 - i]) {
            case 'A' : rc[i] = 'T';
                       break;
            case 'C' : rc[i] = 'G';
                       break;
            case 'G' : rc[i] = 'C';
                       break;
            case 'T' : rc[i] = 'A';
                       break;
            case 'N' : rc[i] = 'N';
                       break;
            default : PyErr_SetString(PyExc_ValueError, "invalid nt");
                      return NULL;
        }
    }
    py_rc = PyUnicode_FromString(rc);
    PyMem_Del(rc);
    return py_rc;
}


static PyObject *
lowQtoN(PyObject *self, PyObject *args)
{
    // define variables
    PyObject *py_newr;
    int minq;
    const char *r, *q;
    size_t rlen, i;

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
            "Same as `dms_tools2.utils.buildReadConsensus` but "
            "`r2` should be reverse-complemented prior to call."},
    {"alignSubamplicon", alignSubamplicon, METH_VARARGS,
            "Same as `dms_tools2.utils.alignSubamplicon`."},
    {"lowQtoN", lowQtoN, METH_VARARGS,
            "Same as `dms_tools2.utils.lowQtoN`."},
    {"reverseComplement", reverseComplement, METH_VARARGS,
            "Same as `dms_tools2.utils.reverseComplement`."},
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

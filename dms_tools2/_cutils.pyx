"""
====================
_cutils
====================

`cython` implementations of some functions in `utils`.

This module implements faster `cython` versions of some of
the functions in the `utils` module of ``dms_tools2``.

The recommended way to call these functions is **not** directly,
but rather via their pure Python implementation in `utils` using
the `use_cutils` argument.
"""


import numpy
cimport numpy
import dms_tools2


# define char * such that _NTS[i] == dms_tools2.NTS[i]
cdef bytes _NTS_bytes = ''.join(dms_tools2.NTS).encode('UTF-8')
cdef char* _NTS = _NTS_bytes


def buildReadConsensus(list reads, int minreads, float minconcur):
    """Faster implementation of `utils.buildReadConsensus`.

    See the doc string of `utils.buildReadConsensus` for
    details on what this function does.
    """
    cdef:
        int i, j
        int maxlen = max(map(len, reads))
        int n_nts = len(dms_tools2.NTS)
        str r
        bytes r_bytes
        char *r_char
        numpy.ndarray[numpy.int_t, ndim=2] counts = numpy.zeros(
                [maxlen, n_nts], dtype=numpy.int)

    for r in reads:
        r_bytes = r.encode('UTF-8')
        r_char = r_bytes
        for i in range(len(r)):
            if r_char[i] != 'N':
                for j in range(n_nts):
                    if r_char[i] == _NTS[j]:
                        counts[i][j] += 1
                        break
                else:
                    raise ValueError("Invalid nucleotide {0}".format(
                            chr(r_char[i])))
    cdef:
        numpy.ndarray[numpy.int_t, ndim=1] ntot = counts.sum(axis=1)
        numpy.ndarray[numpy.int_t, ndim=1] xmax = counts.argmax(axis=1)
    consensus = numpy.full(shape=maxlen, fill_value=chr(_NTS[0]))
    for i in range(1, n_nts):
        numpy.place(consensus, xmax == i, chr(_NTS[i]))
    with numpy.errstate(divide='ignore', invalid='ignore'):
        consensus[((ntot < minreads) |  ((counts.max(axis=1).astype('float') /
                ntot) < minconcur))] = 'N'
    return ''.join(consensus)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

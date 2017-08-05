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


# define str such that _NTS[i] == dms_tools2.NTS[i]
cdef str _NTS = ''.join(dms_tools2.NTS)


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
        cdef numpy.ndarray[numpy.int_t, ndim=2] counts = numpy.zeros(
                [maxlen, n_nts], dtype=numpy.int)

    for r in reads:
        for i in range(len(r)):
            if r[i] != 'N':
                for j in range(n_nts):
                    if r[i] == _NTS[j]:
                        counts[i][j] += 1
                        break
                else:
                    raise ValueError("Invalid nucleotide {0}".format(r[i]))
    cdef:
        numpy.ndarray[numpy.int_t, ndim=1] ntot = counts.sum(axis=1)
        numpy.ndarray[numpy.int_t, ndim=1] xmax = counts.argmax(axis=1)
    consensus = numpy.full(shape=maxlen, fill_value=_NTS[0])
    for i in range(1, n_nts):
        numpy.place(consensus, xmax == i, _NTS[i])
    with numpy.errstate(divide='ignore', invalid='ignore'):
        consensus[((ntot < minreads) |  ((counts.max(axis=1).astype('float') /
                ntot) < minconcur))] = 'N'
    return ''.join(consensus)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

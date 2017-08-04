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


def buildReadConsensus(list reads, int minreads, float minconcur):
    """Faster implementation of `utils.buildReadConsensus`.

    See the doc string of `utils.buildReadConsensus` for
    details on what this function does.
    """
    cdef:
        int i, maxlen, ntot, lenr, nmax

    readlens = list(map(len, reads))
    maxlen = max(readlens)
    consensus = []
    for i in range(maxlen):
        counts = {}
        for (r, lenr) in zip(reads, readlens):
            if lenr > i:
                x = r[i]
                if x != 'N':
                    if x in counts:
                        counts[x] += 1
                    else:
                        counts[x] = 1
        ntot = sum(counts.values())
        if ntot < minreads:
            consensus.append('N')
        else:
            (nmax, xmax) = sorted([(n, x) for (x, n) in counts.items()])[-1]
            if nmax / float(ntot) >= minconcur:
                consensus.append(xmax)
            else:
                consensus.append('N')
    return ''.join(consensus)

if __name__ == '__main__':
    import doctest
    doctest.testmod()

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
        int i, ntot, nmax
        int maxlen = max(map(len, reads))

    counts = [{} for i in range(maxlen)]
    for r in reads:
        for i in range(len(r)):
            x = r[i]
            if x != 'N':
                if x in counts[i]:
                    counts[i][x] += 1
                else:
                    counts[i][x] = 1
    consensus = []
    for i in range(maxlen):
        ntot = sum(counts[i].values())
        if ntot < minreads:
            consensus.append('N')
        else:
            (nmax, xmax) = sorted([(n, x) for (x, n) in counts[i].items()])[-1]
            if nmax / float(ntot) >= minconcur:
                consensus.append(xmax)
            else:
                consensus.append('N')
    return ''.join(consensus)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

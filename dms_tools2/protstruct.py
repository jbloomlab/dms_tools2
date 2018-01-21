"""
===========
protstruct
===========
Module for operations related to protein structures.
"""

import os
import itertools
import tempfile
import collections

import natsort
import numpy
import Bio.PDB


def atomDist(a1, a2):
    """Calculates distance between two atoms.

    Args:
        `a1`, `a2` (`Bio.PDB.Atom.Atom` objects)
            The two atoms.

    Returns:
        The distance between the two atoms.

    >>> Atom = Bio.PDB.Atom.Atom
    >>> a1 = Atom(name='CA', coord=(1.0, 2.0, 3.0), bfactor=1, occupancy=1,
    ...           altloc='', fullname='CA', serial_number=1, element='C')
    >>> a2 = Atom(name='CA', coord=(1.5, 3.0, 4.5), bfactor=1, occupancy=1,
    ...           altloc='', fullname='CA', serial_number=2, element='C')
    >>> round(atomDist(a1, a2), 3)
    1.871
    """
    diff_vector = numpy.array(a1.coord) - numpy.array(a2.coord)
    return numpy.sqrt(numpy.sum(diff_vector * diff_vector))


def residueDist(r1, r2, atom=None):
    """Calculates distance between closest atoms in two residues.

    Args:
        `r1`, `r2` (`Bio.PDB.Residue.Residue` objects)
            The two residues.
        `atom` (`None` or str)
            If `None`, computes closest distance between any atoms.
            If a string, only computes closest distance among atoms
            of this type (e.g., `CA` for alpha carbon distances).

    Returns:
        Distance between closest atoms in the two residues.

    >>> Residue = Bio.PDB.Residue.Residue
    >>> Atom = Bio.PDB.Atom.Atom
    >>> r1 = Residue('r1', 'r1', 'A')
    >>> a1 = Atom(name='CA', coord=(1.0, 2.0, 2.5), bfactor=1, occupancy=1,
    ...           altloc='', fullname='CA', serial_number=1, element='C')
    >>> a2 = Atom(name='N', coord=(1.5, 3.0, 4.5), bfactor=1, occupancy=1,
    ...           altloc='', fullname='N', serial_number=2, element='N')
    >>> r1.add(a1)
    >>> r1.add(a2)
    >>> r2 = Residue('r2', 'r2', 'A')
    >>> a3 = Atom(name='CA', coord=(2.0, 4.0, 6.0), bfactor=1, occupancy=1,
    ...           altloc='', fullname='CA', serial_number=1, element='C')
    >>> r2.add(a3)
    >>> round(residueDist(r1, r2), 3) == 1.871
    True
    >>> round(residueDist(r1, r2, atom='CA'), 3) == 4.153
    True
    """
    if atom is None:
        d = min([atomDist(a1, a2) for (a1, a2) in itertools.product(r1, r2)])
    else:
        d = atomDist(r1[atom], r2[atom])
    return d


def distMatrix(pdbfile, chains, dist_type, equivchains={}, ignore_hetero=True):
    """Residue-residue distance matrix for protein or homo-oligomer.

    Args:
        `pdbfile` (str)
            Name of PDB file.
        `chains` (str or list of str)
            Chain in `pdbfile` for which we compute residue distances,
            or list of chains if the residues span several chains 
            (as for proteins like HA which are cleaved to subunits).
        `dist_type` (str)
            Distances to measure. Use `CA` for alpha carbon distances,
            and `any` for nearest distances of any atom in residues.
        `equivchains` (dict)
            If the structure is a homo-oligomer, each chain may
            have other equivalent chains. In this case, for each
            `chain`, `equivchains[chain]` should be list of equivalent
            chains that we also include in the distance calculations.
        `ignorehetero` (bool)
            Ignore hetero-residues, and only consider protein ones.    

    Returns:
        The 2-tuple `(residues, distmatrix)` where `residues` is a
        list of the residue names (as strings) and `distmatrix` is a 
        `numpy.ndarray` with element `distmatrix[i, j]` giving
        the distance between residue `residues[i]` and `residues[j]`.

    Here is an example computation of distances between two residues
    spanning two chains:

    >>> with tempfile.NamedTemporaryFile(mode='w') as f:
    ...     n = f.write(
    ...      '\\n'.join([
    ...      'ATOM   4633  N   VAL X 505A     57.621  44.297  43.089  1.00 96.43           N',
    ...      'ATOM   4634  CA  VAL X 505A     58.278  45.594  43.147  1.00 84.49           C',
    ...      'ATOM   4635  C   VAL X 505A     59.779  45.434  42.943  1.00 82.23           C',
    ...      'ATOM   4636  O   VAL X 505A     60.339  44.372  43.218  1.00 91.07           O',
    ...      'ATOM   4637  CB  VAL X 505A     58.010  46.304  44.493  1.00 98.34           C',
    ...      'ATOM   4638  CG1 VAL X 505A     58.732  47.644  44.544  1.00 94.96           C',
    ...      'ATOM   4639  CG2 VAL X 505A     56.511  46.483  44.713  1.00 92.78           C',
    ...      'ATOM      1  N   VAL A 518      46.814  16.139  29.171  1.00 92.74           N',
    ...      'ATOM      2  CA  VAL A 518      46.514  15.640  27.833  1.00 99.45           C',
    ...      'ATOM      3  C   VAL A 518      47.047  16.605  26.764  1.00102.76           C',
    ...      'ATOM      4  O   VAL A 518      47.281  17.785  27.034  1.00 86.02           O',
    ...      'ATOM      5  CB  VAL A 518      44.993  15.424  27.640  1.00 95.35           C',
    ...      'ATOM      6  CG1 VAL A 518      44.729  14.428  26.515  1.00 70.29           C',
    ...      'ATOM      7  CG2 VAL A 518      44.357  14.931  28.935  1.00 97.65           C',
    ...      'ATOM      8  CA  VAL Y 505A     47.514  16.640  28.833  1.00 99.45           C',
    ...      ]))
    ...     f.flush()
    ...     (residues, ca_dist) = distMatrix(f.name, ['A', 'X'], 'CA')
    ...     (residues, any_dist) = distMatrix(f.name, ['A', 'X'], 'any')
    ...     (residues, equiv_dist) = distMatrix(f.name, ['A', 'X'], 'CA', {'X':['Y']})
    >>> residues
    ['505A', '518']
    >>> numpy.allclose(ca_dist, numpy.array(
    ...    [[0, 35.639], [35.639, 0]]), atol=1e-3)
    True
    >>> numpy.allclose(any_dist, numpy.array(
    ...    [[0, 32.674], [32.674, 0]]), atol=1e-3)
    True
    >>> numpy.allclose(equiv_dist, numpy.array(
    ...    [[0, 1.732], [1.732, 0]]), atol=1e-3)
    True
    """
    if isinstance(chains, str):
        chains = list(chains)

    # read PDB and take first model (typically there is only one unless NMR)
    structure = Bio.PDB.PDBParser().get_structure('id', pdbfile)
    model = structure[0]

    # lists residue objects 
    res_objs = []
    # keyed by residue objects, lists equivalents in other chains
    equiv_res_objs = collections.defaultdict(list) 
    for chain in chains:
        for res in model[chain]:
            res_id = res.get_id()
            if res_id[0].strip() and ignore_hetero:
                continue
            res_objs.append(res)
            if equivchains and chain in equivchains:
                for otherchain in equivchains[chain]:
                    equiv_res_objs[res] += [res2 for res2 in
                            model[otherchain] if res2.get_id() == res_id]

    # get residue names and sort 
    residues = ['{0}{1}'.format(n, i).strip() for (het, n, i) in
            map(lambda x: x.get_id(), res_objs)]
    assert len(residues) == len(set(residues)), "duplicate residue names"
    sortindex = natsort.index_natsorted(residues)
    residues = natsort.order_by_index(residues, sortindex)
    res_objs = natsort.order_by_index(res_objs, sortindex)

    # build up distance matrix
    if dist_type == 'any':
        atom = None
    elif dist_type == 'CA':
        atom = 'CA'
    else:
        raise ValueError("invalid dist_type: {0}".format(dist_type))
    distmatrix = numpy.zeros((len(residues), len(residues)))
    for ((i, ri), (j, rj)) in itertools.combinations(enumerate(res_objs), 2):
        dist = min(
                [residueDist(ri, rj, atom=atom)] +
                [residueDist(ri, rk, atom=atom) for rk in equiv_res_objs[rj]] +
                [residueDist(rk, rj, atom=atom) for rk in equiv_res_objs[ri]]
                )
        distmatrix[i, j] = distmatrix[j, i] = dist

    return (residues, distmatrix)



if __name__ == '__main__':
    import doctest
    doctest.testmod()

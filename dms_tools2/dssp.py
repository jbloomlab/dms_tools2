"""
===========
dssp
===========
Process output from `dssp <http://swift.cmbi.ru.nl/gv/dssp/>`_.

`dssp <http://swift.cmbi.ru.nl/gv/dssp/>`_ can be used to calculate
secondary structure and solvent accessibility information from a
PDB structure. This module can process that output.
"""

import os
import re
import pandas
import Bio.PDB


#: max accessible surface area (square angstroms) for amino acids, from
#: `Tien et al (2013) <https://doi.org/10.1371/journal.pone.0080635>`_.
#: `MAX_ASA_TIEN[a]` is the max surface area for amino acid `a`.
MAX_ASA_TIEN = {'A':129.0, 'R':274.0, 'N':195.0, 'D':193.0, 'C':167.0,
                'E':223.0, 'Q':225.0, 'G':104.0, 'H':224.0, 'I':197.0,
                'L':201.0, 'K':236.0, 'M':224.0, 'F':240.0, 'P':159.0,
                'S':155.0, 'T':172.0, 'W':285.0, 'Y':263.0, 'V':174.0}


def processDSSP(dsspfile, chain=None, max_asa=MAX_ASA_TIEN):
    """Get secondary structure and solvent accessibility from ``dssp``.

    `dssp <http://swift.cmbi.ru.nl/gv/dssp/>`_ is a program
    that calculates secondary structure and absolute solvent
    accessibility from a PDB file.

    This function processes the text output provided by the
    `dssp webserver <http://swift.cmbi.ru.nl/gv/dssp/>`_, at
    least given the format of that output as of Sept-4-2017.

    It returns a `pandas.DataFrame` that gives the secondary
    structure and solvent accessibility for each residue in the
    ``dssp`` output.

    Args:
        `dsspfile` (str)
            Name of text file containing ``dssp`` output.

        `chain` (str or `None`)
            If the PDB file analyzed by ``dssp`` to create `dsspfile`
            has more than one chain, specify the letter code for one
            of those chains with this argument.

        `max_asa` (dict)
            Max surface area for each amino acid in square angstroms.

    Returns:
        A `pandas.DataFrame` with the following columns:

            - `site`: residue number for all sites in `dsspfile`.

            - `amino_acid`: amino acid identity of site in `dsspfile`.

            - `ASA`: absolute solvent accessibility of the residue.

            - `RSA`: relative solvent accessibility of the residue.

            - `SS`: ``dssp`` secondary structure code, one of:

                - *G*: 3-10 helix

                - *H*: alpha helix

                - *I*: pi helix

                - *B*: beta bridge

                - *E*: beta bulge

                - *T*: turn

                - *S*: high curvature

                - *-*: loop

            - `SS_class`: broader secondary structure class:

                - *helix*: `SS` value of *G*, *H*, or *I*

                - *strand*: `SS` value of *B* or *E*

                - *loop*: any of the other `SS` values.
    """
    dssp_cys = re.compile('[a-z]')
    d_dssp = Bio.PDB.make_dssp_dict(dsspfile)[0]
    chains = set([chainid for (chainid, r) in d_dssp.keys()])
    if chain is None:
        assert len(chains) == 1, "chain is None, but multiple chains"
        chain = list(chains)[0]
    elif chain not in chains:
        raise ValueError("Invalid chain {0}".format(chain))
    d_df = {'site':[],
            'amino_acid':[],
            'ASA':[],
            'RSA':[],
            'SS':[],
            'SS_class':[],
            }
    for ((chainid, r), tup) in d_dssp.items():
        if chainid == chain:
            (tmp_aa, ss, acc) = tup[ : 3]
            if dssp_cys.match(tmp_aa):
                aa = 'C'
            else:
                aa = tmp_aa
            if r[2] and not r[2].isspace():
                # site has letter suffix
                d_df['site'].append(str(r[1]) + r[2].strip())
            else:
                d_df['site'].append(r[1])
            d_df['amino_acid'].append(aa)
            d_df['ASA'].append(acc)
            d_df['RSA'].append(acc / float(max_asa[aa]))
            d_df['SS'].append(ss)
            if ss in ['G', 'H', 'I']:
                d_df['SS_class'].append('helix')
            elif ss in ['B', 'E']:
                d_df['SS_class'].append('strand')
            elif ss in ['T', 'S', '-']:
                d_df['SS_class'].append('loop')
            else:
                raise ValueError("invalid SS of {0}".format(ss))
    return pandas.DataFrame(d_df)



if __name__ == '__main__':
    import doctest
    doctest.testmod()

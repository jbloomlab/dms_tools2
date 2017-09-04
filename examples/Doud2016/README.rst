===============================================================
Doud 2016 deep mutational scanning of influenza hemagglutinin
===============================================================

This subdirectory contains a re-analysis using `dms_tools2 <https://jbloomlab.github.io/dms_tools2/>`_ of the deep mutational scanning of the A/WSN/1933 (H1N1) hemagglutinin by `Doud and Bloom (2016) <http://www.mdpi.com/1999-4915/8/6/155>`_.

The analysis is performed by the Jupyter notebook `analysis_notebook.ipynb <analysis_notebook.ipynb>`_.

The requisite input data are in the ``./data/`` subdirectory:

    * `./data/WSN-HA.fasta <./data/WSN-HA.fasta>`_ contains the wildtype WSN HA coding sequence.

    * `./data/originalDoud2016prefs <./data/originalDoud2016prefs>`_ contains the original amino-acid preferences estimated by `Doud and Bloom (2016) <http://www.mdpi.com/1999-4915/8/6/155>`_ using the older `dms_tools <https://jbloomlab.github.io/dms_tools/>`_ software.

    * `./data/1RVX_trimer_sequentialnumbering.pdb <./data/1RVX_trimer_sequentialnumbering.pdb>`_ is the H1 HA trimer in `PDB 1RVX <http://www.rcsb.org/pdb/explore.do?structureId=1rvx>`_, re-numbered with `PDB Goodies <http://dicsoft2.physics.iisc.ernet.in/pdbgoodies/inputpage.html>`_ to a sequential numbering scheme starting with 1 at the N-terminal Met.

    * `./data/1RVX_trimer_sequentialnumbering.dssp <./data/1RVX_trimer_sequentialnumbering.dssp>`_ is the downloaded result of running the `DSSP webserver <http://swift.cmbi.ru.nl/gv/dssp/>`_ on `./data/1RVX_trimer_sequentialnumbering.pdb <./data/1RVX_trimer_sequentialnumbering.pdb>`_, and so holds absolute solvent accessibility and secondary structure information.

    * `./data/HA_alignment.fasta/ <./data/HA_alignment.fasta>`_ is an alignment of H1 HA sequences. It is **not** the same alignment used by `Doud and Bloom (2016) <http://www.mdpi.com/1999-4915/8/6/155>`_; it is substantially smaller for one thing.

====================================================================
Dingens 2017 mutational antigenic profiling of HIV Env
====================================================================

This subdirectory contains a re-analysis using `dms_tools2 <https://jbloomlab.github.io/dms_tools2/>`_ of the mutational antigenic profiling of the Env protein from the BF520 strain of HIV by `Dingens, Haddox, Overbaugh, and Bloom (2017) <https://doi.org/10.1016/j.chom.2017.05.003>`_.

The analysis is performed by the Jupyter notebook `analysis_notebook.ipynb <analysis_notebook.ipynb>`_.

The requisite input data are in the ``./data/`` subdirectory:

    * `./data/BF520c2-Env.fasta <./data/BF520c2-Env.fasta>`_ contains the wildtype coding sequence of the BF520 Env used in the experiment.

    * `./data/BF520c2_to_HXB2.csv <./data/BF520c2_to_HXB2.csv>`_ gives the mapping from sequential 1, 2, ... numbering of the BF520 protein sequence to the standard `HXB2 numbering scheme <https://www.hiv.lanl.gov/content/sequence/HIV/REVIEWS/HXB2.html>`_ for HIV Env.

    * `./data/LAIalaninescandata_logbase2.csv <.data/LAIalaninescandata_logbase2.csv>`_  and `./data/JRCSFalaninescandata_logbase2_115foldlimit.csv <JRCSFalaninescandata_logbase2_115foldlimit.csv>`_ give the results of the partial mostly alanine scanning mutagenesis of the LAI and JRCSF strains by `Falkowska et al, 2014 <http://www.sciencedirect.com/science/article/pii/S107476131400123X>`_. These fiels give the log2 fold change in IC50 relative to wildtype, using the lower maximum endpoint as the maximum effect in both files.

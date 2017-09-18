====================================================================
Doud 2017 mutational antigenic profiling of influenza hemagglutinin
====================================================================

This subdirectory contains a re-analysis using `dms_tools2 <https://jbloomlab.github.io/dms_tools2/>`_ of the mutational antigenic profiling of the A/WSN/1933 (H1N1) hemagglutinin by `Doud, Hensley, and Bloom (2017) <http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006271>`_.

The analysis is performed by the Jupyter notebook `analysis_notebook.ipynb <analysis_notebook.ipynb>`_.

The requisite input data are in the ``./data/`` subdirectory:

    * `./data/WSN-HA.fasta <./data/WSN-HA.fasta>`_ contains the wildtype WSN HA coding sequence.

    * `./data/known_escape.csv <./data/known_escape.csv>`_ is a file that contains the sites of the classically identified escape mutants from the four antibodies used in this study. These sites are numbered in sequential 1, 2, ... numbering of HA. It is the same set compiled in `Doud, Hensley, and Bloom (2017) <http://journals.plos.org/plospathogens/article?id=10.1371/journal.ppat.1006271>`_. The set was compiled from the following references:

        - escape mutants in `Caton et al, 1982 <https://www.ncbi.nlm.nih.gov/pubmed/6186384>`_

        - escape mutants in `Das et al, 2013 <https://www.ncbi.nlm.nih.gov/pubmed/23498956>`_

        - antibody naming scheme in `Magadan et al, 2013 <http://jvi.asm.org/content/87/17/9742.full>`_

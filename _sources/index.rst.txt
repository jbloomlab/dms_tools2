.. dms_tools2 documentation master file, created by
   sphinx-quickstart on Fri Jul 28 17:39:58 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for ``dms_tools2``
======================================

`dms_tools2`_ is a software package for analyzing `deep mutational scanning <https://www.ncbi.nlm.nih.gov/pubmed/25075907>`_ data.
It is tailored to analyze libraries created using `comprehensive codon mutagenesis`_ of protein-coding genes, and perform analyses that are common to the `Bloom lab`_, such as:

    * Process deep-sequencing data from :ref:`bcsubamp`.

    * Estimate the :ref:`prefs` of a protein.

    * Estimate :ref:`diffsel` when imposing a pressure such as antibody selection.

    * Estimate the :ref:`fracsurvive` for each mutation when imposing a strong pressure such as antibody selection.

    * Processing PacBio CCS sequencing of defined amplicons using the :py:mod:`dms_tools2.pacbio` module.

The `dms_tools2 source code`_ is freely available under a `GPLv3`_ license.
`dms_tools2`_ is a re-write of the earlier `dms_tools <http://github.com/jbloomlab/dms_tools>`_ package.

If you use `dms_tools2`_ for your work, **please cite the references** in :ref:`citations`.

:ref:`installation` of `dms_tools2`_ installs a collection of :ref:`programs` as well as a :ref:`api`.

The easiest way to learn to use `dms_tools2`_ is to look at the :ref:`examples`.

Contents
----------
.. toctree::
   :maxdepth: 2

   installation
   examples
   bcsubamp
   prefs
   diffsel
   fracsurvive
   programs
   api
   citations


Indices and tables
---------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. include:: weblinks.txt

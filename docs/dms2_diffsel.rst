.. _dms2_diffsel:

==========================================
``dms2_diffsel``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_diffsel`` program processes files giving the number of observed counts of characters in a selected and mock-selected condition to estimate :ref:`diffsel`.

If you have multiple related replicates or samples (or even if you hvae just one), you should probably use the :ref:`dms2_batch_diffsel` program rather than running ``dms2_diffsel`` directly.
This is because :ref:`dms2_batch_diffsel` runs ``dms2_diffsel``, but then also makes some nice summary plots.

.. _prefs_commandlineusage:

Command-line usage
----------------------------------------
.. argparse::
   :module: dms_tools2.parseargs
   :func: diffselParser
   :prog: dms2_diffsel

   \-\-sel
    The counts files have the format of the files created by programs such as ``dms2_bcsubamp``. 
    Specifically, they must have the following columns: 'site', 'wildtype', and then a column for each possible character (e.g., codon).

   \-\-name
    The `Output files`_ will have a prefix equal to the name specified here.
    This name should only contain letters, numbers, and dashes.

   \-\-indir
    This option can be useful if the counts files are found in a common directory. 
    Instead of repeatedly listing that directory name, you can just provide it here.

.. _diffsel_outputfiles:

Output files
--------------
The output files all have the prefix specified by ``--outdir`` and ``--name``.
For instance, if you use ``--outdir results --name replicate-1``, then the output files will have the prefix ``./results/replicate-1`` and the suffixes described below.

Here are the specific output files:

Log file
+++++++++++
This file has the suffix ``.log``. 
It is a text file that logs the progress of the program.

.. include:: weblinks.txt

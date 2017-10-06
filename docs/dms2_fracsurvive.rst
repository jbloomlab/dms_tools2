.. _dms2_fracsurvive:

==========================================
``dms2_fracsurvive``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_fracsurvive`` program processes files giving the number of observed counts of characters in a selected and mock-selected condition along with a measurement of the overall fraction of the library surviving the selection to estimate the :ref:`fracsurvive` for each mutation.

If you have multiple related replicates or samples (or even if you have just one), you should probably use the :ref:`dms2_batch_fracsurvive` program rather than running ``dms2_fracsurvive`` directly.
This is because :ref:`dms2_batch_fracsurvive` runs ``dms2_fracsurvive``, but then also makes some nice summary plots.

.. _fracsurvive_commandlineusage:

Command-line usage
----------------------------------------
.. argparse::
   :module: dms_tools2.parseargs
   :func: fracsurviveParser
   :prog: dms2_fracsurvive

   \-\-sel
    The counts files have the format of the files created by programs such as ``dms2_bcsubamp``. 
    Specifically, they must have the following columns: 'site', 'wildtype', and then a column for each possible character (e.g., codon).

   \-\-name
    The `Output files`_ will have a prefix equal to the name specified here.
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

   \-\-indir
    This option can be useful if the counts files are found in a common directory. 
    Instead of repeatedly listing that directory name, you can just provide it here.

.. _fracsurvive_outputfiles:

Output files
--------------
The output files all have the prefix specified by ``--outdir`` and ``--name``.
For instance, if you use ``--outdir results --name replicate-1``, then the output files will have the prefix ``./results/replicate-1`` and the suffixes described below.

Here are the specific output files:

Log file
+++++++++++
This file has the suffix ``.log``. 
It is a text file that logs the progress of the program.

Mutation fraction surviving file
+++++++++++++++++++++++++++++++++++++++
This file has the suffix ``_mutfracsurvive.csv``.
It gives the fraction surviving for each mutation at each site, which is the :math:`F_{r,x}` value defined in Equation :eq:`fracsurvive` of the :ref:`fracsurvive` section.
Note that the quantity is calculated for the wildtype as well as the mutant characters at each site.
Note also that if you are using ``--aboveavg yes`` then these are the fraction surviving **above the library average**, denoted as :math:`F_{r,x}^{\rm{aboveavg}}` in Equation :eq:`fracsurviveaboveavg` of the :ref:`fracsurvive` section.
If ``--mincounts`` is greater than zero, the fraction surviving may be undefined for some mutations due to low counts, and any such undefined values are also shown as `NaN`.

Here are the first and last few lines of a ``_mutfracsurvive.csv`` file::

    site,wildtype,mutation,mutfracsurvive
    156,G,S,0.8189280293912643
    146,N,D,0.626080490632122
    157,K,S,0.5933429890043687
    158,S,A,0.5723610875357631
    ...
    540,L,C,0.0016521655105078295
    175,P,G,0.0013556328151850907
    545,S,T,0.001310686952133144
    490,E,D,0.001016678753248357

Note that the file is sorted from largest to smallest fraction surviving.

Site fraction surviving file
+++++++++++++++++++++++++++++++++
This file has the suffix ``_sitefracsurvive.csv``
It gives several measures that summarize the fraction surviving each site.
All values in the ``_sitefracsurvive.csv`` file can be calculated from the values in the ``_mutfracsurvive.csv`` file, but the program outputs both files to make things simpler for the user.

Specifically, it gives the following quantities:

* `avgfracsurvive` is the average of the mutation fraction surviving values. If any of the mutation fraction surviving values are `NaN` (which can happen if you use ``--mincounts``), they are **not** included in this average.

* `maxfracsurvive` is the **maximum** mutation fraction surviving taken over all **non**-wildtype characters for each site.

Here are the first and last lines of a ``_sitefracsurvive.csv`` file::

    site,avgfracsurvive,maxfracsurvive
    153,0.2841412228950011,0.54440652255242
    157,0.25281575693873276,0.5933429890043687
    136,0.16447395413487326,0.33671209426557835
    156,0.13827268994538547,0.8189280293912643
    ...
    210,0.013137991003787609,0.022298531591000946
    170,0.011505865316217256,0.029289469175867944
    176,0.010814303017871948,0.02157678717020969
    175,0.008361156202363286,0.027496593984557578

If all mutations at a site have a mutation fraction surviving of `NaN` (which can be the case if ``--mincounts`` is > 0), then the site values are reported as `NaN`.

.. include:: weblinks.txt

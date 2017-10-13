.. _dms2_batch_fracsurvive:

==========================================
``dms2_batch_fracsurvive``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_batch_fracsurvive`` program can be used to estimate :ref:`fracsurvive` for each mutation following a stringent selection such as antbody treatment.

The ``dms2_batch_fracsurvive`` program runs :ref:`dms2_fracsurvive` for each sample listed in a batch file specified by ``--batchfile``.
Specifically, as described in :ref:`batch_diffsel_commandlineusage`, you can specify a few sample-specific arguments in the ``--batchfile``.
All other arguments are specified using the normal option syntax (e.g., ``--indir INDIR``) and are shared between all samples specified in ``--batchfile``.
The result is the output for each individual run of :ref:`dms2_fracsurvive` plus the summary plots described in `Output files`_.
It then creates the summary plots described in `Output files`_.

Because ``dms2_batch_fracsurvive`` simply runs :ref:`dms2_fracsurvive` on each sample specified by the ``--batchfile`` argument described below, see the :ref:`dms2_fracsurvive` :ref:`fracsurvive_commandlineusage` for details that are helpful for understanding some of the arguments in the ``dms2_batch_fracsurvive`` :ref:`batch_fracsurvive_commandlineusage` below.

.. _batch_fracsurvive_commandlineusage:

Command-line usage
---------------------------------------------

.. argparse::
   :module: dms_tools2.parseargs
   :func: batch_fracsurviveParser
   :prog: dms2_batch_fracsurvive

   \-\-batchfile
    Each of the arguments `name`, `sel`, `mock`, `libfracsurvive`, and optionally `err` gives the value of the same parameter passed to ``dms2_fracsurvive``. 
    If `group` is being used, then the `group` is pre-pended to `name` for that sample.
    In addition, `group` is used to organize output for similar runs that should be grouped when calculating means / medians and plotting.

    If you are running with no error-control counts, then do **not** specify ``--err``.

   \-\-summaryprefix
    As detailed in `Output files` below, ``dms2_batch_fracsurvive`` creates a variety of plots summarizing the output.
    These files are in the directory specified by ``--outdir``, and have the prefix specified here.
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

Output files
--------------
Running ``dms2_batch_fracsurvive`` produces output files in the directory specified by ``--outdir``.

Results for each sample
++++++++++++++++++++++++++
The program :ref:`dms2_fracsurvive` is run on each sample specified by ``--batchfile``, so you will create all of the :ref:`dms2_fracsurvive` :ref:`fracsurvive_outputfiles`.

If you are using the `group` entry in ``--batchfile``, then for each sample we create a name by pre-pending the group to the name. 
For instance, if ``--batchfile`` is::

    group,name,sel,mock
    antibody-1,replicate-1,sel_1_1,mock_1
    antibody-1,replicate-2,sel_1_2,mock_2
    antibody-2,replicate-1,sel_2_1,mock_1
    antibody-2,replicate-2,sel_2_2,mock_2

then the output files for the individual samples will have prefixes like ``antibody-1-replicate-1_*``, ``antibody-1-replicate-2_*``, etc.

On the other hand, if ``--batchfile`` does not specify groups, then the name for each sample is just given by the `name` column. 
So if ``--batchfile`` is::

    name,sel,mock
    replicate-1,sel_1_1,mock_1
    replicate-2,sel_1_2,mock_2

then the output files will have prefixes like ``replicate-1_*``, ``replicate-2_*``.

Mean and median fraction surviving
++++++++++++++++++++++++++++++++++++++++
The program computes the mean and median fraction surviving values for each group (if there are groups), or for all samples.
Note that the means and medians are computed on the mutation fraction surviving values, and then the site values are computed from these mean / median mutation selections.
The files are in the same format as those created by :ref:`dms2_fracsurvive`.

For instance, for the first example ``--batchfile`` in the section above (the one with a `group` column), we would get the following files if we used ``--summaryprefix summary``::

    summary_antibody-1-meanmutfracsurvive.csv
    summary_antibody-2-meanmutfracsurvive.csv
    summary_antibody-1-medianmutfracsurvive.csv
    summary_antibody-2-medianmutfracsurvive.csv
    summary_antibody-1-meansitefracsurvive.csv
    summary_antibody-2-meansitefracsurvive.csv
    summary_antibody-1-mediansitefracsurvive.csv
    summary_antibody-2-mediansitefracsurvive.csv

For the second example ``--batchfile`` (the one without a `group` column), we would get the following files::

    summary_meanmutfracsurvive.csv
    summary_medianmutfracsurvive.csv
    summary_meansitefracsurvive.csv
    summary_mediansitefracsurvive.csv

It is often useful to visualize the mean or median `mutfracsurvive` files with :ref:`dms2_logoplot`.

Correlation plots
+++++++++++++++++++++++++++++++++
Scatter plots are created that show the correlations among samples within the same `group`, or among all samples if there are not any groups.

Separate plots are made for the `mutfracsurvive`, the `avgfracsurvive` (averaged across all mutations at each site), and the maximum `mutfracsurvive` at each site.
The names will have the form::

    summary_antibody-1-mutfracsurvivecorr.pdf
    summary_antibody-1-avgfracsurvivecorr.pdf
    summary_antibody-1-maxfracsurvivecorr.pdf

Fracsurvive plots
+++++++++++++++++++
Plots are made that show the site average and maximum fraction surviving as a function of the primary sequence. 
These plots show the mean and median values for each group, and are faceted by `group` (if there are groups).
If you run with ``--summaryprefix summary``, then the plots will be:

    - `avgfracsurvive`: files ``summary_meanavgfracsurvive.pdf`` and ``summary_medianavgfracsurvive.pdf`` show the average `mutfracsurvive` across all mutations for each site.

    - `maxfracsurvive`: files ``summary_meanmaxfracsurvive.pdf`` and ``summary_medianmaxfracsurvive.pdf`` show the maximum `mutfracsurvive` across all mutations for each site.

Log file
++++++++++++
A log file is created that summarizes the output.
For instance, if you run ``dms2_batch_fracsurvive`` with the arguments ``--outdir results --summaryprefix summary`` then the log will be ``./results/summary.log``.


.. include:: weblinks.txt

.. _dms2_batch_diffsel:

==========================================
``dms2_batch_diffsel``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_batch_diffsel`` program can be used to estimate :ref:`diffsel`.

The ``dms2_batch_diffsel`` program simply runs :ref:`dms2_diffsel` for each sample listed in a batch file specified by ``--batchfile``.
Specifically, as described in :ref:`batch_diffsel_commandlineusage`, you can specify a few sample-specific arguments in the ``--batchfile``.
All other arguments are specified using the normal option syntax (e.g., ``--indir INDIR``) and are shared between all samples specified in ``--batchfile``.
The result is the output for each individual run of :ref:`dms2_diffsel` plus the summary plots described in `Output files`_.
It then creates the summary plots described in `Output files`_.

The `Doud2017 example`_ to illustrates the usage of ``dms2_batch_diffsel`` on a real dataset.

Because ``dms2_batch_diffsel`` simply runs ``dms2_diffsel`` on each sample specfied by the ``--batchfile`` argument described below, see the ``dms2_diffsel`` :ref:`diffsel_commandlineusage` for details that are helpful for understanding some of the arguments in the ``dms2_batch_diffsel`` :ref:`batch_diffsel_commandlineusage` below.

.. _batch_diffsel_commandlineusage:

Command-line usage
---------------------------------------------

.. argparse::
   :module: dms_tools2.parseargs
   :func: batch_diffselParser
   :prog: dms2_batch_diffsel

   \-\-batchfile
    Each of the arguments `name`, `sel`, `mock`, and optionally `err` gives the value of the same parameter passed to ``dms2_diffsel``. 
    If `group` is being used, then the `group` is pre-pended to `name` for that sample.
    In addition, `group` is used to organize output for similar runs that should be grouped when calculating means / medians and plotting.

    If you are running with no error-control counts, then do **not** specify ``--err``.

   \-\-summaryprefix
    As detailed in `Output files` below, ``dms2_batch_diffsel`` creates a variety of plots summarizing the output.
    These files are in the directory specified by ``--outdir``, and have the prefix specified here.
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

Output files
--------------
Running ``dms2_batch_diffsel`` produces output files in the directory specified by ``--outdir``.

Results for each sample
++++++++++++++++++++++++++
The program ``dms2_diffsel`` is run on each sample specified by ``--batchfile``, so you will create all of the ``dms2_diffsel`` :ref:`diffsel_outputfiles`.

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

Mean and median differential selection
++++++++++++++++++++++++++++++++++++++++
The program computes the mean and median differential selection for each group (if there are groups), or for all samples.
Note that the means and medians are computed on the mutation differential selection, and then the site differential selection values are computed from these mean / median mutation differential selections.
The files are in the same format as those created by :ref:`dms2_diffsel`.

For instance, for the first example ``--batchfile`` in the section above (the one with a `group` column), we would get the following files if we used ``--summaryprefix summary``::

    summary_antibody-1-meanmutdiffsel.csv
    summary_antibody-2-meanmutdiffsel.csv
    summary_antibody-1-medianmutdiffsel.csv
    summary_antibody-2-medianmutdiffsel.csv
    summary_antibody-1-meansitediffsel.csv
    summary_antibody-2-meansitediffsel.csv
    summary_antibody-1-mediansitediffsel.csv

For the second example ``--batchfile`` (the one without a `group` column), we would get the following files::

    summary_meanmutdiffsel.csv
    summary_medianmutdiffsel.csv
    summary_meansitediffsel.csv
    summary_mediansitediffsel.csv

It is often useful to visualize the mean or median mutdiffsel files with :ref:`dms2_logoplot`.

Correlation plots
+++++++++++++++++++++++++++++++++
Scatter plots are created that show the correlations among samples within the same `group`, or among all samples if there are not any groups.

Separate plots are made for the mutdiffsel, the absolute sitediffsel, the positive sitediffsel, and the maximum mutdiffsel at each site.
The names will have the form::

    summary_antibody-1-mutdiffselcorr.pdf
    summary_antibody-1-absolutesitediffselcorr.pdf
    summary_antibody-1-positivesitediffselcorr.pdf
    summary_antibody-1-maxmutdiffselcorr.pdf

For examples of these plots, see the `Doud2017 example`_.

Diffsel plots
+++++++++++++++
Plots are made that show the differential selection as a function of the primary sequence. 
These plots show the mean and median values for each group, and are faceted by `group` (if there are groups).
If you run with ``--summaryprefix summary``, then the plots will be:

    - total sitediffsel: files ``summary_meantotaldiffsel.pdf`` and ``summary_mediantotaldiffsel.pdf`` show both positive and negative sitediffsel.

    - positive sitediffsel: files ``summary_meanpositivediffsel.pdf`` and ``summary_medianpositivediffsel.pdf`` show just positive sitediffsel.

    - minmax sitediffsel: files ``summary_meanminmaxdiffsel.pdf`` and ``summary_medianminmaxdiffsel.pdf`` show minimum and maximum mutdiffsel for each site.

    - max sitediffsel: files ``summary_meanmaxdiffsel.pdf`` and ``summary_medianmaxdiffsel.pdf`` show maximum mutdiffsel for each site.

For examples of these plots, see the `Doud2017 example`_.

Log file
++++++++++++
A log file is created that summarizes the output.
For instance, if you run ``dms2_batch_diffsel`` with the arguments ``--outdir results --summaryprefix summary`` then the log will be ``./results/summary.log``.


.. include:: weblinks.txt

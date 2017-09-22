.. _dms2_batch_prefs:

==========================================
``dms2_batch_prefs``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_batch_prefs`` program processes files giving the number of observed counts of characters pre- and post-selection to estimate :ref:`prefs`.

The ``dms2_batch_prefs`` program simply runs :ref:`dms2_prefs` for each sample listed in a batch file specified by ``--batchfile``.
Specifically, as described in :ref:`batch_prefs_commandlineusage`, you can specify a few sample-specific arguments in the ``--batchfile``.
All other arguments are specified using the normal option syntax (e.g., ``--indir INDIR``) and are shared between all samples specified in ``--batchfile``.
The result is the output for each individual run of :ref:`dms2_prefs` plus the summary plots described in `Output files`_.

The ``dms2_batch_prefs`` program simply runs :ref:`dms2_prefs` for each sample listed in a batch file.
It then creates the summary plots described in `Output files`_.

The `Doud2016 example`_ to illustrates the usage of ``dms2_batch_prefs`` on a real dataset.

Because ``dms2_batch_prefs`` simply runs ``dms2_prefs`` on each sample specfied by the ``--batchfile`` argument described below, see the ``dms2_prefs`` :ref:`prefs_commandlineusage` for details that are helpful for understanding many of the arguments in the ``dms2_batch_prefs`` :ref:`batch_prefs_commandlineusage` below.

.. _batch_prefs_commandlineusage:

Command-line usage
---------------------------------------------

.. argparse::
   :module: dms_tools2.parseargs
   :func: batch_prefsParser
   :prog: dms2_batch_prefs

   \-\-batchfile
    Each of the arguments `name`, `pre`, and `post` gives the value of the same parameter as passed to ``dms2_prefs``.

    If you are running with no error-control counts, then do **not** specify ``--err`` or ``--errpre`` / ``--errpost``.

    If you are running with a single error-control for both the pre- and post-selection counts, then specify this counts file with ``--err``.

    If you have different controls for the pre- and post-selection counts, specify them separately with ``--errpre`` and ``--errpost``.

   \-\-summaryprefix
    As detailed in `Output files` below, ``dms2_batch_prefs`` creates a variety of plots summarizing the output.
    These files are in the directory specified by ``--outdir``, and have the prefix specified here.
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

Output files
--------------
Running ``dms2_batch_prefs`` produces output files in the directory specified by ``--outdir``.

Results for each sample
++++++++++++++++++++++++++
The program ``dms2_prefs`` is run on each sample specified by ``--batchfile``, so you will create all of the ``dms2_prefs`` :ref:`prefs_outputfiles`.

Average preferences
+++++++++++++++++++++
A file is created that holds the preferences averaged across all samples in ``--batchfile``.
This file has the prefix specified by ``--summaryprefix``. 
For instance, if you run ``dms2_batch_prefs`` with the arguments ``--outdir results --summaryprefix summary`` then the plot will be ``./results/summary_avgprefs.csv``.
It has the same format as the preferences files created by ``dms2_prefs``.

Correlation plot
+++++++++++++++++++++++++++++++++
A plot is created that summarizes the correlation between the preferences for each sample in ``--batchfile``.
This plot has the prefix specified by ``--summaryprefix``. 
For instance, if you run ``dms2_batch_prefs`` with the arguments ``--outdir results --summaryprefix summary`` then the plot will be ``./results/summary_prefscorr.pdf``.
An example of this plot is in the `Doud2016 example`_.

Log file
++++++++++++
A log file is created that summarizes the output.
For instance, if you run ``dms2_batch_prefs`` with the arguments ``--outdir results --summaryprefix summary`` then the log will be ``./results/summary.log``.

Program run time
---------------------------
As described in the :ref:`prefs_runtime` section for ``dms2_prefs``, each iteration of that program can take a while to run.
So obviously running it multiple times with ``dms2_batch_prefs`` will take even longer.
The time can be reduced by specifying more CPUs to use with ``--ncpus``.

.. include:: weblinks.txt

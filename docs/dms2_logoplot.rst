.. _dms2_logoplot:

==========================================
``dms2_logoplot``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_logoplot`` program makes logo-plot visualizations of the data.
It uses a slightly modified version of `weblogo`_ to make the logo plots themselves.

If you run it using ``--prefs`` to specify the input file, then the plot will visualize :ref:`prefs`.
An example of the type of plot created for :ref:`prefs` is shown in the `Doud2016 example`_.

If you run it using ``--diffsel`` to specify the input file, then the plot will visualize :ref:`diffsel`.
An example of the type of plot created for :ref:`diffsel` is shown in the `Doud2017 example`_.

.. _batch_logoplot_commandlineusage:

Command-line usage
---------------------------------------------

.. argparse::
   :module: dms_tools2.parseargs
   :func: logoplotParser
   :prog: dms2_logoplot


Output files
--------------
Running ``dms2_logoplot`` produces output files in the directory specified by ``--outdir``, and with the prefix specified by ``--name``.

There will be a log file with the suffix ``.log`` summarizing the program's progress.

If you run with ``--prefs``, then the logo plot will be in a file with the suffix ``_prefs.pdf``.
An example of such a logo plot is in the `Doud2016 example`_.

If you run with ``--diffsel``, then the logo plot will be a file with the suffix ``_diffsel.pdf``.
An example of such a logo plot is in the `Doud2017 example`_.


.. include:: weblinks.txt

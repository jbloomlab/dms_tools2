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

If you run it using ``--muteffects`` to specify the input file, then the plot will visualize the logarithm base 2 of the mutational effect.
You can calculate mutational effects from :ref:`prefs` using the function :func:`dms_tools2.prefs.prefsToMutFromWtEffects`.
The mutational effect calculated by this function is just the ratio of the preference for the mutant amino acid over the preference for the wildtype amino acid.

If you run it using ``--diffsel`` to specify the input file, then the plot will visualize :ref:`diffsel`.

If you run it using ``--fracsurvive`` to specify the input file, then the plot will visualize the :ref:`fracsurvive` for each mutation.

If you run it using ``--diffprefs`` to specify the input file, then the plot will show the difference between preferences, showing negative and positive values.

.. _logoplot_commandlineusage:

Command-line usage
---------------------------------------------

.. argparse::
   :module: dms_tools2.parseargs
   :func: logoplotParser
   :prog: dms2_logoplot

   \-\-name
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

Output files
--------------
Running ``dms2_logoplot`` produces output files in the directory specified by ``--outdir``, and with the prefix specified by ``--name``.

There will be a log file with the suffix ``.log`` summarizing the program's progress.

If you run with ``--prefs``, then the logo plot will be in a file with the suffix ``_prefs.pdf``.
An example of such a logo plot is in the `Doud2016 example`_.

If you run with ``--diffsel``, then the logo plot will be a file with the suffix ``_diffsel.pdf``.
An example of such a logo plot is in the `Doud2017 example`_.

If you run with ``--fracsurvive``, then the logo plot will be a file with the suffix ``_fracsurvive.pdf``.

If you run with ``--muteffects``, then the logo plot will be a file with the suffix ``_muteffects.pdf``.

If you run with ``--diffprefs``, then the logo plot will be a file with the suffix ``_diffprefs.pdf``.


.. include:: weblinks.txt

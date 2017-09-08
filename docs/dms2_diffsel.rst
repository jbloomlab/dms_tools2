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

.. _diffsel_commandlineusage:

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
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

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

Mutation differential selection file
+++++++++++++++++++++++++++++++++++++++
This file has the suffix ``_mutdiffsel.csv``.
It gives the differential selection for each mutation at each site, which is the :math:`s_{r,x}` value defined in Equation :eq:`mutdiffsel` of the :ref:`diffsel` section.
The mutation differential selection values are shown as ``NaN`` for the wildtype identity at a site.
If ``--mincounts`` is greater than zero, the differential selection may also be undefined for some mutations due to low counts, and any such undefined differential selection values are also shown as ``NaN``.

Here are the first and last few lines of a ``_mutdiffsel.csv`` file::

    site,wildtype,mutation,mutdiffsel
    156,G,S,8.20616611420281
    157,K,S,7.843970369736835
    146,N,D,7.839003488125009
    157,K,I,7.636912452846413
    153,S,I,7.618894947835398
    ...
    560,Q,Q,NaN
    561,C,C,NaN
    562,R,R,NaN
    563,I,I,NaN
    564,C,C,NaN
    565,I,I,NaN

Note that the file is sorted from largest to smallest mutation differential selection, with ``NaN`` values last.

Site differential selection file
+++++++++++++++++++++++++++++++++
This file has the suffix ``_sitediffsel.csv``
It gives several measures that summarize the differential selection at each site.
All values in the ``_sitediffsel.csv`` file can be calculated from the values in the ``_mutdiffsel.csv`` file, but we output both files to make things simpler for the user.

Specifically, it gives the following quantities:

1. ``abs_diffsel`` is the total mutation differential selection (both positive and negative) at a site, as defined in Equation :eq:`abs_diffsel` in the :ref:`diffsel` section.

2. ``positive_diffsel`` is the total positive mutation differential selection at a site, as defined in Equation :eq:`positive_diffsel` in the :ref:`diffsel` section.

3. ``negative_diffsel`` is the total negative mutation differential selection at a site, as defined in Equation :eq:`negative_diffsel` in the :ref:`diffsel` section.

4. ``max_diffsel`` is the maximum mutation differential selection at a site, as defined in Equation :eq:`max_diffsel` in the :ref:`diffsel` section.

5. ``min_diffsel`` is the minimum mutation differential selection at a site, as defined in Equation :eq:`min_diffsel` in the :ref:`diffsel` section.

Here are the first and last lines of a ``_sitediffsel.csv`` file::

    site,abs_diffsel,positive_diffsel,negative_diffsel,max_diffsel,min_diffsel
    157,103.49320489879341,103.49320489879341,0.0,7.843970369736835,0.0
    153,83.14675113695142,81.29729853121766,-1.8494526057337632,7.618894947835398,-1.5180701665199172
    156,50.6281426671937,50.6281426671937,0.0,8.20616611420281,0.0
    158,44.54232019822232,44.09637071052413,-0.4459494876981922,6.721775476391176,-0.3980554768143438
    175,32.69746741931557,0.0,-32.69746741931557,0.0,-3.5655146998011067
    ...
    274,16.78563070533967,3.9335896820646017,-12.85204102327507,2.5127254463901822,-1.5263088393038542
    171,16.687556395725117,5.582191979487671,-11.105364416237446,5.286815935040484,-1.3671039380789383
    229,16.589450103457924,0.0,-16.589450103457924,0.0,-2.4947190310966607
    299,16.3418067952726,0.0,-16.3418067952726,0.0,-1.596540281615842
    499,16.33889299748386,6.3537025445869935,-9.985190452896864,3.564608664594661,-2.218441347118039

If all mutations at a site have a mutation differential selection of ``NaN`` (which can be the case if ``--mincounts`` is > 0), then the site differential selection values are reported as 0.

.. include:: weblinks.txt

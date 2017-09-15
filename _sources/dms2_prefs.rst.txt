.. _dms2_prefs:

==========================================
``dms2_prefs``
==========================================

.. contents::
   :local:

Overview
-------------
The ``dms2_prefs`` program processes files giving the number of observed counts of characters pre- and post-selection to estimate :ref:`prefs`.

If you have multiple replicates, you should probably use the :ref:`dms2_batch_prefs` program rather than running ``dms2_prefs`` directly.

.. _prefs_commandlineusage:

Command-line usage
----------------------------------------
.. argparse::
   :module: dms_tools2.parseargs
   :func: prefsParser
   :prog: dms2_prefs

   \-\-pre
    The counts files have the format of the files created by programs such as ``dms2_bcsubamp``. 
    Specifically, they must have the following columns: 'site', 'wildtype', and then a column for each possible character (e.g., codon).

   \-\-name
    The `Output files`_ will have a prefix equal to the name specified here.
    This name should only contain letters, numbers, dashes, and spaces.
    Underscores are **not** allowed as they are a LaTex special character.

   \-\-indir
    This option can be useful if the counts files are found in a common directory. 
    Instead of repeatedly listing that directory name, you can just provide it here.

.. _prefs_outputfiles:

Output files
--------------
The output files all have the prefix specified by ``--outdir`` and ``--name``.
For instance, if you use ``--outdir results --name replicate-1``, then the output files will have the prefix ``./results/replicate-1`` and the suffixes described below.

Here are the specific output files:

Log file
+++++++++++
This file has the suffix ``.log``. 
It is a text file that logs the progress of the program.

Preferences file
++++++++++++++++++++++
This file has the suffix ``_prefs.csv``. 
It gives the estimate preference for each character at each site. 
For instance::

    site,A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y
    1,0.05201542494870646,0.05006892743588129,0.04935929298955449,0.04561969825792428,0.05027584699870549,0.049473267870167204,0.04966234700207467,0.05035563145797863,0.05180290264027603,0.052502385368045426,0.05323583603885511,0.046134643822643394,0.051011235226750974,0.052905549301723524,0.04290757945900044,0.049609243914675485,0.06315116938502258,0.04129711732370402,0.051893451567002306,0.04671844899130815
    2,0.009456820014302609,0.0887864717361851,0.03202933078899705,0.01243744054076179,0.012598634010240075,0.016079011606326327,0.14360828098089565,0.05270848041464168,0.05479041395821132,0.059805747035781634,0.10546956692273086,0.03356178777566189,0.012037714482087113,0.08503842285654405,0.017824696236264426,0.027336538959474643,0.05028936751898547,0.029519124422321352,0.01629464236061828,0.1403275073789687
    3,0.094394784492365,0.033233499951948485,0.10037681454416572,0.041772952245424946,0.01871075286571138,0.010914843906391419,0.01994461441568695,0.09430640509845868,0.010261045290749045,0.050955385392754314,0.06764316761334091,0.06593302352530313,0.047625012474641924,0.017370598629944167,0.1082951339123566,0.04003184839931041,0.07144380858649375,0.026212403552398438,0.02646517359744569,0.05410873150510903
    4,0.07817657215908004,0.03148741643399614,0.005538443259083886,0.018851757050952038,0.0034453072574090094,0.030655060310952557,0.03370373802129379,0.023488641120853936,0.05342118049856918,0.05175840113766944,0.2235830210977376,0.07104192962903758,0.03487046604114975,0.0796424680240337,0.052235719104467615,0.02309884775188897,0.05227025898510587,0.04266732483424344,0.04636513033841905,0.04369831694405645

.. _prefs_runtime:

Program run time
---------------------------
If you run ``dms_prefs`` with ``--method ratio`` then it will run very quickly.

If you run it with ``--method bayesian`` then the runtime will be somewhat longer due to the MCMC.
Exactly how long depends on whether you are using error controls for the counts (the ``--err`` option).
If you use different files for the pre- and post-selection error controls, and are using ``--chartype codon_to_aa`` then the program will typically take about 4 or 5 hours if you give it 4 CPUs.
If you give it more CPUs, or using the same (or no) error control for pre- and post-selection, then it will be faster.

.. include:: weblinks.txt

.. _examples:

==========================================
Examples
==========================================

.. contents::
   :local:

Example analyses
------------------
Each analysis is in a `Jupyter notebook`_ in its own subdirectory `on GitHub <https://github.com/jbloomlab/dms_tools2/tree/master/examples>`_.
Specifically:

Deep mutational scanning of influenza hemagglutinin by Doud and Bloom (2016)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
`Doud and Bloom (2016)`_ performed deep mutational scanning of influenza hemagglutinin using :ref:`bcsubamp` to obtain high sequencing accuracy.
You can see a `Jupyter notebook`_ that analyzes their data by clicking here: `Doud2016 example`_.

Mutational antigenic profiling of influenza hemagglutinin by Doud et al (2017)
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
`Doud et al (2017)`_ performed mutational antigenic profiling of influenza hemagglutinin against a set of strain-specific antibodies, using :ref:`bcsubamp` to obtain high sequencing accuracy.
You can see a `Jupyter notebook`_ that analyzes their data by clicking here: `Doud2017 example`_.

Downloading and running the example notebooks
-------------------------------------------------
Each example analysis is available on the GitHub repository that hosts the `dms_tools2 source code`_. 
Each analysis is in the form of a `Jupyter notebook`_, and other necessary input data is also provided.
`Navigate here <https://github.com/jbloomlab/dms_tools2/tree/master/examples>`_ to access these notebooks and the associated data.

To run an example `Jupyter notebook`_, you will need to install ``jupyter`` for Python 3.
If you don't have that already installed, you can install it via::

    pip install jupyter --user

assuming that you have already set up your paths as described for the :ref:`installation` of `dms_tools2`_.
Note that as for :ref:`installation` of `dms_tools2`_, if ``pip`` points to the Python 2 version on your computer, then you need to use ``pip3``.
Most of the notebooks also have other requirements (such as `fastq-dump`_ to download the deep sequencing data from the `Sequence Read Archive`_), and these dependencies are described in the notebooks.

To run a notebook interactively, type::

    jupyter notebook analysis_notebook.ipynb

To run a notebook in command-line mode, run::

    jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 analysis_notebook.ipynb

You should **not** open the notebook while it is being run in command-line mode.

If you are in the `Bloom lab`_ and using the Hutch server, you can run a notebook via `slurm`_ by creating a file ``run.sbatch`` that looks like this::

    #!/bin/sh
    #SBATCH
    #PBS -l walltime=24:00:00
    jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=-1 analysis_notebook.ipynb

and then running the file with::

    sbatch -p largenode -c 14 --mem=140000 run.sbatch

where the number of CPUs requested with ``-c`` should match the number specified for multi-processor operations in the notebook.


.. include:: weblinks.txt

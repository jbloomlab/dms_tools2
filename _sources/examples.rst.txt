.. _examples:

==========================================
Examples
==========================================

.. contents::
   :local:

Example analyses
------------------
Each example analysis is a published study by the `Bloom lab`_.
Each analysis is in its own subdirectory `in the examples folder on GitHub <https://github.com/jbloomlab/dms_tools2/tree/master/examples>`_.
Typically there is a `Jupyter notebook`_ that performs the analysis.

Here are the analyses:

Doud and Bloom (2016)
+++++++++++++++++++++++
 - **Summary**: Deep mutational scanning of influenza hemagglutinin using :ref:`bcsubamp`, then estimation of the :ref:`prefs`.
 - **Citation**: `Doud and Bloom (2016)`_
 - **Jupyter notebook**: `Doud2016 example`_

Doud et al (2017)
++++++++++++++++++++++
 - **Summary**: Mutational antigenic profiling of influenza hemagglutinin against strain-specific antibodies using :ref:`bcsubamp`, then estimation of the :ref:`diffsel` from each antibody
 - **Citation**: `Doud et al (2017)`_
 - **Jupyter notebook**:`Doud2017 example`_.

Dingens et al (2017)
++++++++++++++++++++++
 - **Summary**: Mutational antigenic profiling of HIV Env (BF520 strain) against PGT151 antibody using :ref:`bcsubamp`, and then estimation of the :ref:`diffsel`.
 - **Citation**: `Dingens et al (2017)`_
 - **Jupyter notebook**: `Dingens2017 example`_.

Doud et al (2018)
+++++++++++++++++++
 - **Summary**: Mutational antigenic profiling of influenza hemagglutinin against broadly neutralizing and strain-specific antibodies using :ref:`bcsubamp`, then estimation of the :ref:`fracsurvive`.
 - **Citation**: `Doud et al (2018)`_
 - **Jupyter notebook**: `Doud2018 example`_.

Haddox et al (2018)
+++++++++++++++++++++
 - **Summary**: Deep mutational scanning of two HIV Envs using :ref:`bcsubamp`, then comparison of the :ref:`prefs`.
 - **Citation**: `Haddox et al (2018)`_
 - **Jupyter notebook**: `Haddox2018 example`_.

Dingens et al (2018)
+++++++++++++++++++++
 - **Summary**: Mutational antigenic profiling of HIV Env against fusion-peptide antibodies usring :ref:`bcsubamp`, then estimation of the :ref:`diffsel`.
 - **Citations**: `Dingens et al (2018)`_
 - **Jupyter notebook**: `Dingens2018 example`_

Lee et al (2018)
+++++++++++++++++
 - **Summary**: Deep mutational scanning of H3 influenza hemagglutinin using :ref:`bcsubamp`, then estimation of the :ref:`prefs`.
 - **Citations**: `Lee et al (2018)`_
 - **Jupyter notebook**: `Lee2018 example`_

Running the examples
-------------------------------------------------
`Navigate here <https://github.com/jbloomlab/dms_tools2/tree/master/examples>`_ to access the example analyses and require input data.

To run a `Jupyter notebook`_, you need to install ``jupyter``.
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

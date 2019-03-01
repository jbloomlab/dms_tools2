.. _installation:

================================
Installation
================================

.. contents::
   :local:

Make sure you have Python **3.6** or higher
--------------------------------------------
`dms_tools2`_ requires `Python`_ **3.6** or higher.
Most computers have both Python 2 and Python 3 installed.
You can check your current version of `Python`_ (and the Python installation program `pip`_) with::

    python --version

and::

    pip --version

If these commands indicate that your current default versions are Python 2, then you need to either re-set the defaults to be Python 3 or use the Python 3 specific commands on your computer (which are likely to be ``python3`` rather than ``python``, and ``pip3`` rather than ``pip``).


Installing with ``pip`` 
----------------------------------------------------

Quick installation
++++++++++++++++++++++++++
If your system already has the appropriate version of ``pip`` and appropriate paths, you can install `dms_tools2`_ with the simple command::

    pip install dms_tools2 --user

If this command fails, then read the instructions below.


Where to install
+++++++++++++++++
You need to figure out where you want to install `dms_tools2`_.
Global installation using ``sudo`` `is not recommended for Python packages in general <http://stackoverflow.com/questions/21055859/what-are-the-risks-of-running-sudo-pip/21056000#21056000>`_.

The simplest solution is to install locally via the ``--user`` option to ``pip``, which by default on Linux will install into the ``~/.local/`` directory.

In order for locally installed programs to be accessible, you need to add ``~/.local/bin/`` to the ``PATH`` variable, and ``~/.local/lib/`` to the ``PYTHONPATH`` variable. If you are using the `bash shell`_, you would do this by adding the following lines to your ``~/.bashrc`` file::

    PATH=$HOME/.local/bin/:$PATH
    export PYTHONPATH=$HOME/.local/lib/python3.6:$PATH

You then want to make sure that your ``~/.bash_profile`` file simple sources your ``~/.bashrc`` file as `described here <http://www.joshstaiger.org/archives/2005/07/bash_profile_vs.html>`_ by making ``~/.bash_profile`` consist of the following contents::

    if [ -f ~/.bashrc ]; then
        source ~/.bashrc
    fi

On Mac OS X, the default directory for ``--user`` may be ``$HOME/Library/Python/x.y/`` rather than ``~/.local/`` where ``x.y`` indicates the version number (e.g., ``3.6``.

Make sure ``pip`` is installed
++++++++++++++++++++++++++++++++++++++++++

Check if you already have `pip`_ installed. You can do this by typing at the command line::

    pip --version

If this command indicates that you have `pip`_ for Python 3.6 or higher, then you can move to the next step. 
If you do not have `pip`_, then you need to install it by following the `instructions here <https://pip.pypa.io/en/stable/installing/>`_.

Use ``pip`` to install ``dms_tools2``
++++++++++++++++++++++++++++++++++++++++++
Once `pip`_ is installed, you can do a local installation with::

    pip install dms_tools2 --user

Using a virtual environment
++++++++++++++++++++++++++++++
The other good option rather than ``--user`` is to use ``pip`` to install into a virtual environment `as described here <http://docs.python-guide.org/en/latest/dev/virtualenvs/>`_

Installing the ``rplot`` module
++++++++++++++++++++++++++++++++++
**Note**: we now recommend that you use the `dmslogo`_ package rather than the `rplot module`_.

As described in the :ref:`api`, there are some useful features
in the `rplot module`_. In order for your installation to
support this module, you need to install a recent version of 
`R <https://www.r-project.org/about.html>`_ and then run your installation
with::

    pip install dms_tools2[rplot] --user


Upgrading with ``pip``
--------------------------------------------------
If you have previously installed `dms_tools2`_ but are not sure that you have the latest version, you can upgrade using `pip`_. To do this for a local installation, use::

    pip install dms_tools2 --user --upgrade


Install from source code
-----------------------------------------------------------------------
You can also install the latest version of the `dms_tools2 source code`_ from GitHub. 

To install from source, first clone the `dms_tools2 source code`_ from GitHub::

    git clone https://github.com/jbloomlab/dms_tools2

Then install locally with::

    cd dms_tools2
    pip install -e . --user

If you have already cloned the repository, you can update the source by::

    cd dms_tools2
    git pull
    pip install -e . --user

If you want to install the `rplot module`_ from source, the command is::

    pip install -e .[rplot] --user

.. _license:

License
-----------
`dms_tools2 source code`_ is available on GitHub under an open-source `GPLv3`_ license. Part of the code utilized by `dms_tools2`_ is based on `weblogo`_, which is licensed under the GPL-compatible BSD 3-clause license.


.. include:: weblinks.txt

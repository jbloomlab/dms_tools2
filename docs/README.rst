===========================
Documentation
===========================

This subdirectory contains the `reStructuredText`_ documentation that can be built with `sphinx`_.

Building the documentation
-----------------------------

To build the documentation, you need to install:

    * `sphinx`_ 
    
    * `sphinx-argparse`_ 

Then simply type::

    make html

and the HTML documentation will be built in ``./_build/html/``.

Notes
--------

Note that the configuration automatically created by ``sphinx-quickstart`` has been modified in the following ways:

    * ``conf.py`` has been modified to read the version and package information from ``../dms_tools2/_metadata.py``

    * ``Makefile`` has been modified to automatically run `sphinx-apidoc`_.


.. _`reStructuredText`: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _`sphinx`: http://sphinx-doc.org/
.. _`sphinx-argparse`: http://sphinx-argparse.readthedocs.org

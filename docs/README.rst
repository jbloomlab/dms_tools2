===========================
Documentation
===========================

This subdirectory contains the `reStructuredText`_ documentation that can be built with `sphinx`_.

Building the documentation
-----------------------------

To build the documentation, you need to install:

    * `sphinx`_ 
    
    * `sphinx-argparse`_ 

    * `nb2plots`_

Then simply type::

    make html

and the HTML documentation will be built in ``./_build/html/``.

Notes on configuration
------------------------

Note that the configuration automatically created by ``sphinx-quickstart`` has been modified in the following ways:

    * ``conf.py`` has been modified to:
    
        - use `sphinx-argparse`_ to enable parsing of command-line program arguments
        
        - read the version and package information from ``../dms_tools2/_metadata.py``

        - specify `numfig = True` to enable figure numbering

    * ``Makefile`` has been modified to automatically run `sphinx-apidoc`_.

Notes on nb2plots
-------------------
`nb2plots`_ was used to generate `codonvariant_analysis_demo.rst <codonvariant_analysis_demo.rst>`_ from a Jupyter notebook.

.. _`reStructuredText`: http://docutils.sourceforge.net/docs/user/rst/quickref.html
.. _`sphinx`: http://sphinx-doc.org/
.. _`sphinx-argparse`: http://sphinx-argparse.readthedocs.org
.. _`sphinx-apidoc`: http://www.sphinx-doc.org/en/stable/man/sphinx-apidoc.html
.. _`nb2plots`: https://matthew-brett.github.io/nb2plots/

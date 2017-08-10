==============
tests
==============

This directory contains some tests for ``dms_tools2``. 

Each ``test_*.py`` file runs a different test, and can be executed with commands like::

    python test_bcsubamplicons.py

The recommended way to run all these tests **and** the doctests in the individual modules is using `pytest`_.
To do this, first install `pytest`_ if it is not already installed.
Then navigate to the **top directory** of ``dms_tools2`` (**not** this ``./tests/`` directory) and run `../test.bash <../test.bash>`_.


.. _`pytest`: https://docs.pytest.org/en/latest/

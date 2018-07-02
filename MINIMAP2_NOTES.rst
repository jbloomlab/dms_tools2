====================================================
Notes on inclusion of `minimap2` as a `git subtree`
====================================================

The `minimap2`_ executable is included as a `git subtree <https://developer.atlassian.com/blog/2015/05/the-power-of-git-subtree/>`_, and then built and installed as package data by `setup.py <setup.py>`_.
This installed version is then used by default by `dms_tools2.minimap2`.
It is installed in a location where it is **not** expected to be in the normal system path to be available at the command line.
This provides a consistently versioned `minimap2`_ within the package.

Version 2.10 of `minimap2`_ was added to ``./minimap2_source/`` via::

    git subtree add --prefix minimap2_source https://github.com/lh3/minimap2 v2.10 --squash

It was then updated to version 2.11 with::

    git subtree pull --prefix minimap2_source https://github.com/lh3/minimap2 v2.11 --squash

All the files in this ``./minimap2_source/`` directory are specified for inclusion in the main package via `MANIFEST.in <MANIFEST.in>`_.

The `setup.py <setup.py>`_ script then takes a somewhat hacky approach to build the ``minimap2`` executable and copy it to ``./dms_tools2/`` under the name ``minimap2_prog`` for inclusion as package data.

.. _`minimap2`: https://github.com/lh3/minimap2

Change Log
===========

2.2.dev0
---------
* added `compareprefs` module

* Added `omega` overlay option to ``dms2_logoplot``

* Fix bug with ``dms2_logoplot`` when using wildtype sequence overlays

* Fix bug with ``--fracsurvivemax 0`` to ``dms2_logoplot``

* Fix bug with handling of disulfide-bonded cysteines in ``dssp`` output.

* Added `colors` option to `plot.plotCorrMatrix`

2.1.0
------
* Added programs and docs for `fracsurvive`.

* Added ``--scalebar`` to ``dms2_logoplot``.

* Add `grouplabel` option and preserve group order for faceted plots by batch programs.

* Handle dependencies without `__version__` attribute

2.0.2
------
* Added ``--sitemask`` option to ``dms2_bcsubamp`` / ``dms2_batch_bcsubamp``.

* Standardized color scheme in ``*_cumulmutcounts.pdf`` plot.

* Ensure naturally sorted average prefs from ``dms2_batch_prefs``.

2.0.1
------
* A few packaging changes for PyPI

2.0.0
--------
This version is a complete re-write of `dms_tools <https://github.com/jbloomlab/dms_tools>`_ version 1.2.2.

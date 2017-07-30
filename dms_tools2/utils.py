"""
===================
utils
===================

Miscellaneous utilities for ``dms_tools2``.
"""


import os
import sys
import time
import platform
import importlib
import dms_tools2


def sessionInfo():
    """Returns string with information about session / packages."""
    s = [
            'Version information:',
            '\tTime and date: {0}'.format(time.asctime()),
            '\tPlatform: {0}'.format(platform.platform()),
            '\tPython version: {0}'.format(
                    sys.version.replace('\n', ' ')),
            '\tdms_tools2 version: {0}'.format(dms_tools2.__version__),
            ]
    for modname in ['Bio']:
        try:
            v = importlib.import_module(modname).__version__
            s.append('\t{0} version: {1}'.format(modname, v))
        except ImportError:
            raise ImportError("Cannot import {0}".format(modname))
    return '\n'.join(s)


def initLogger(logfile, prog, args):
    """Initialize `logging.Logger` for scripts.

    Args:
        `logfile` (str)
            Name of file to which log is written.
        `prog` (str)
            Name of program for which we are logging.
        `args` (dict)
            Program arguments as arg / value pairs.

    Returns:
        An opened and initialized `logging.Logger`.
        Basic information about the program and args
        will have already been written to this logger.
    """
    if os.path.isfile(logfile):
        os.remove(logfile)
    logging.basicConfig(level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s')
    logger = logging.getLogger(prog)
    logfile_handler = logging.FileHandler(logfile)
    logger.addHandler(logfile_handler)
    formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s')
    logfile_handler.setFormatter(formatter)
    logger.info("Beginning execution of {0} in directory {1}\n".format(
            prog, os.getcwd()))
    logger.info("Progress is being logged to {0}".format(logfile))
    logger.info(sessionInfo())
    logger.info("Parsed the following arguments:\n\t{0}".format(
            '\n\t'.join(['{0} = {1}'.format(arg, val) for (arg, val)
            in args.items()])))


if __name__ == '__main__':
    import doctest
    doctest.testmod()

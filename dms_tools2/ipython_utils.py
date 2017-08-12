"""
===================
ipython_utils
===================

Utilities for using ``dms_tools2`` in an `IPython` (`Jupyter`) notebook.
"""


import os
import subprocess
import time
import IPython.display


def showPDF(pdfs, width=None):
    """Displays PDF images side-by-side in `IPython` notebook.

    This function is useful over just directly displaying the PDF
    as many web browsers do not show inline PDFs. This function
    converts the PDF to PNG and displays that. It requires that
    ImageMagick ``convert`` be available for the conversion.

    Args:
        `pdfs` (str or list)
            Filename of a PDF, or a list of such filenames.
            Multiple images are displayed side-by-side.
        `width` (float or int)
            Width of the displayed PDF.
    """
    if not isinstance(pdfs, list):
        pdfs = [pdfs]
    assert all(map(os.path.isfile, pdfs)), "Can not find PDFs:\n{0}".format(
            '\n'.join(pdfs))
    try:
        x = subprocess.check_output(['convert', '--version'], 
                stderr=subprocess.STDOUT)
    except:
        raise RuntimeError("Cannot find 'convert' executable")
    try:
        # get temporary name for PNG file
        png = os.path.join(os.path.dirname(pdfs[0]), '._' +
                os.path.splitext(os.path.basename(pdfs[0]))[0] + '.png')
        x = subprocess.check_output(['convert', '-density', '134', '-trim'] +\
                ['-splice', '50x0'] + pdfs + ['+append', png])
        time.sleep(0.5)
        IPython.display.display(IPython.display.Image(png, width=width))
    finally:
        if os.path.isfile(png):
            os.remove(png)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

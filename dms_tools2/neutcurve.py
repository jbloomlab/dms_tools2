"""
=============
neutcurve
=============
Module for fitting and analyzing neutralization curves.
"""


import scipy
import scipy.optimize


class fourParamLogistic:
    """Fit a 4-parameter logistic neutralization curve.

    Fits :math:`f(c) = b + \\frac{t - b}{1 + (c/m)^s}`
    where :math:`f(c)` is the fraction surviving at
    concentration :math:`c`, :math:`m` is the midpoint
    of the neutralization curve, :math:`t` is the top
    value (e.g., 1), :math:`b` is the bottom value (e.g., 0),
    and :math:`s` is the slope of the curve.

    This documentation is written for the case when 
    :math:`f(c)` is the fraction **surviving**, meaning that
    :math:`f(c)` gets smaller as :math:`c` gets larger.
    This should lead you to fit :math:`s > 0`.
    If you instead use :math:`f(c)` as the fraction neutralized,
    you will fit :math:`s < 1`. In this case, the interpretation
    of `top` and `bottom` will be reversed. However, this class
    has not been extensively tested in that scenario.

    To initialize and fit an object, provide the following:
        `cs` (array-like)
            Concentrations for which we have measurements
            to fit.
        `fs` (array-like)
            Same length as `cs`, with `fs[i]` giving
            fraction surviving at `cs[i]`.
        `fixbottom` (bool)
            If `False`, we fit the bottom value of the
            curve. If you instead want to fix it,
            provide the number to fix it to (typically 0).
            Fix this number if you think that the neutralization
            goes to completion at high enough antibody.
            concentration. For some viruses (such as HIV),
            the neutralization never goes to completion as 
            there are some resistant viruses (such as due
            to glycan heterogeneity) -- if that is the case,
            you want to fit the bottom rather than fix it to 0.
        `fixtop` (`False` or a float)
            If `False`, we fit the top value of the
            curve. If you instead want to fix it,
            provide the number to fix it to (typically 1).
            Usually you do **not** want to fit this value,
            as the fraction surviving should always be one
            at sufficiently low antibody concentration.

    After initialization, you can access the following
    attributes:
        `cs` (numpy array)
            Array of the concentration measurements we have
            fitted, sorted from lowest to highest. 
        `fs` (numpy array)
            Array of the fraction surviving measurements
            corresponding to the entries in `cs`.
        `midpoint` (float)
            Midpoint of curve, :math:`m` in equation above.
        `slope` (float)
            Hill slope of curve, :math:`s` in equation above.
        `bottom` (float)
            Bottom of curve (value as :math:`c` gets large),
            :math:`b` in equation above.
        `top` (float)
            Top of curve (value as :math:`c` get small),
            :math:`t` in equation above.
        `ic50` (float)
            Concentration :math:`c` where :math:`f(c) = 0.5`,
            or `None` if this value is outside the range
            of fitted concentrations in `cs`.

    As an example, we first simulate some data with known
    parameter values:

    >>> m = 0.03
    >>> s = 1.9
    >>> b = 0.1
    >>> t = 1.0
    >>> cs = [0.002 * 2**x for x in range(9)]
    >>> fs = [fourParamLogistic.evaluate(c, m, s, b, t) for c in cs]

    Now we fit to these data, and then confirm that the
    fitted values are close to the ones used for the
    simulation:

    >>> neut = fourParamLogistic(cs, fs)
    >>> scipy.allclose(neut.midpoint, m)
    True
    >>> scipy.allclose(neut.slope, s)
    True
    >>> scipy.allclose(neut.top, t)
    True
    >>> scipy.allclose(neut.bottom, b)
    True

    Note how since we fit the curve to simulated data where
    the bottom was 0.1 rather than 0, the midpoint and IC50
    are different. Specifically, the IC50 is larger than the
    midpoint, as you have to go past the midpoint to get
    down to value of 0.5 fraction surviving.

    >>> neut.ic50 > neut.midpoint
    True
    >>> scipy.allclose(neut.ic50, 0.0337385586)
    True
    >>> scipy.allclose(0.5, neut.fracsurvive(neut.ic50))
    True
    >>> neut.fracsurvive(neut.midpoint) > 0.5
    True

    Now here is an example where we constrain both the top
    and the bottom (to 1 and 0, respectively) and fit
    the curve. Now the midpoint and IC50 are the same:

    >>> b2 = 0
    >>> t2 = 1
    >>> fs2 = [fourParamLogistic.evaluate(c, m, s, b2, t2) for c in cs]
    >>> neut2 = fourParamLogistic(cs, fs2, fixbottom=b2)
    >>> scipy.allclose(neut2.midpoint, m)
    True
    >>> scipy.allclose(neut2.ic50, m)
    True

    Now let's fit to concentrations that are all **less**
    than the midpoint, so that we never get anywhere
    close to complete neutralization. In this case,
    the estimated IC50 is unreliable, and so will be
    returned as `None`:
    
    >>> cs3 = [1e-5 * 2**x for x in range(7)]
    >>> (cs3[-1] < m)
    True
    >>> fs3 = [fourParamLogistic.evaluate(c, m, s, b2, t2) for c in cs3]
    >>> neut2 = fourParamLogistic(cs3, fs3, fixbottom=b2)
    >>> neut2.ic50 is None
    True

    """

    def __init__(self, cs, fs, fixbottom=False, fixtop=1):
        """See main class docstring."""
        # get data into arrays sorted by concentration
        self.cs = scipy.array(cs)
        self.fs = scipy.array(fs)
        self.fs = self.fs[self.cs.argsort()]
        self.cs = self.cs[self.cs.argsort()]

        # make initial guess for slope to have the right sign
        if self.fs[0] >= self.fs[-1]:
            self.slope = 1.5
        else:
            self.slope = -1.5

        # make initial guess for top and bottom
        if fixtop is False:
            if self.slope > 0:
                self.top = self.fs.max()
            else:
                self.top = self.fs.min()
        else:
            assert isinstance(fixtop, (int, float))
            self.top = fixtop
        if fixbottom is False:
            if self.slope > 0:
                self.bottom = self.fs.min()
            else:
                self.bottom = self.fs.max()
        else:
            assert isinstance(fixbottom, (int, float))
            self.bottom = fixbottom

        # make initial guess for midpoint
        midval = (self.top - self.bottom) / 2.0
        if (self.fs > midval).all():
            if self.slope > 0:
                self.midpoint = self.cs[-1]
            else:
                self.midpoint = self.cs[0]
        elif (self.fs <= midval).all():
            if self.slope > 0:
                self.midpoint = self.cs[0]
            else:
                self.midpoint = self.cs[-1]
        else:
            # get first index where f crosses midpoint
            i = scipy.argmax((self.fs > midval)[:-1] !=
                    (self.fs > midval)[1:])
            assert (self.fs[i] > midval) != (self.fs[i + 1] > midval)
            self.midpoint = (self.cs[i] + self.cs[i + 1]) / 2.0

        # set up function and initial guesses
        if fixtop is False and fixbottom is False:
            initguess = [self.midpoint, self.slope, self.bottom, self.top]
            func = self.evaluate
        elif fixtop is False:
            initguess = [self.midpoint, self.slope, self.top]
            def func(c, m, s, t):
                return self.evaluate(c, m, s, self.bottom, t)
        elif fixbottom is False:
            initguess = [self.midpoint, self.slope, self.bottom]
            def func(c, m, s, b):
                return self.evaluate(c, m, s, b, self.top)
        else:
            initguess = [self.midpoint, self.slope]
            def func(c, m, s):
                return self.evaluate(c, m, s, self.bottom, self.top)

        (popt, pcov) = scipy.optimize.curve_fit(
                func,
                self.cs,
                self.fs,
                initguess
                )

        if fixtop is False and fixbottom is False:
            (self.midpoint, self.slope, self.top, self.bottom) = popt
        elif fixtop is False:
            (self.midpoint, self.slope, self.top) = popt
        elif fixbottom is False:
            (self.midpoint, self.slope, self.bottom) = popt
        else:
            (self.midpoint, self.slope) = popt


    @property
    def ic50(self):
        """IC50 value.
        
        Concentration where `fracsurvive` is 0.5.
        Equals `midpoint` if `top = 1` and `bottom = 0`.

        Is `None` if IC50 is outside range of concentrations
        `cs` used for the fitting.

        Calculated from:
        :math:`0.5 = b + \\frac{t - b}{1 + (ic50/m)^s}`,
        which solves to 
        :math:`ic50 = m \\times \left(\\frac{t - 0.5}{0.5 - b}\\right)^{1/s}`
        """
        if (self.fs > 0.5).all() or (self.fs < 0.5).all():
            return None
        ic50 = (self.midpoint * ((self.top - 0.5) / 
                (0.5 - self.bottom))**(1.0 / self.slope))
        if ic50 > self.cs[-1] or ic50 < self.cs[0]:
            return None
        else:
            return ic50


    def fracsurvive(self, c):
        """Fraction surviving at `c` for fitted parameters."""
        return self.evaluate(c, self.midpoint, self.slope,
                self.bottom, self.top)


    @staticmethod
    def evaluate(c, m, s, b, t):
        """Returns :math:`f(c) = b + \\frac{t - b}{1 + (c/m)^s}`."""
        return b + (t - b) / (1 + (c / m)**s)



if __name__ == '__main__':
    import doctest
    doctest.testmod()

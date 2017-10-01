.. _fracsurvive:

==========================================
Fraction surviving
==========================================

.. contents::
   :local:
   :depth: 2

Overview
-------------
In some cases, you might subject a deep mutational scanning library to a harsh
selection that purges most mutations.
An example of such an experiment is selecting a virus library with a monoclonal
antibody -- most variants will be inactivated by the antibody, but those with mutations
in the epitope will survive.
If the experiment measures the fraction of the entire library
that survives the selection as well as using sequencing to identify the frequencies
of each mutation among the surviving variants, then it is possible to compute the *fraction of variants with each mutation that survive the selection.*
This quantity is expected to range from zero (if all variants with a specific mutation are inactivated) to one (if all variants with a specific mutation survive the selection).

The `dms_tools2`_ software contains programs to estimate the fraction of each mutation surviving from the overall fraction of the library that survives the selection and the deep sequencing data.
Specifically, you can make these estimates using :ref:`dms2_batch_fracsurvive`, which in turn calls :ref:`dms2_fracsurvive`.
These programs create visualizations of the results, and you can also visualize the fraction surviving with :ref:`dms2_logoplot`.

Definition of fraction surviving
---------------------------------------------
We will use :math:`F_{r,x}` to denote the fraction of variants with character :math:`x` (e.g., an amino acid) at site :math:`r` that survive the selection.

In order to calculate :math:`F_{r,x}`, we need three pieces of information:

    1. The **overall** fraction of the entire library that survives the selection. We will denote this quantity as :math:`\gamma`, and it should range from zero to one. For instance, :math:`\gamma` might be calculated by using qPCR to determine the number of virions that infect cells in the presence of antibody versus the number that infect cells in the absence of antibody -- the ratio of these two numbers would be the overall fraction of the entire library the survives the selection. You could also compute :math:`\gamma` directly from deep-sequencing if your library contains some internal control variant that you know will not be affected by the selection (e.g., a variant with a totally orthogonal envelope protein that is not recognized by the antibody).

    2. The number of deep sequencing counts of each character :math:`x` at each site :math:`r` in both the selected (e.g., antibody neutralized) and mock-selected samples. Denote these quantities as :math:`n_{r,x}^{\rm{selected}}` and :math:`n_{r,x}^{\rm{mock}}`, respectively. Also define variables that give the total number of counts at each site for each condition as :math:`N_r^{\rm{selected}} = \sum_x n_{r,x}^{\rm{selected}}` and :math:`N_r^{\rm{mock}} = \sum_x n_{r,x}^{\rm{mock}}`.

Given these pieces of information, the estimated fraction of variants in the selected and mock-selected libraries that contain character :math:`x` at site :math:`r` are calculated as:

.. math::
   :label: fracsurvive_rho

   \rho_{r,x}^{\rm{selected}} &=& \frac{n_{r,x}^{\rm{selected}} + f_{r, \rm{selected}} \times P}{N_r^{\rm{selected}} + f_{r, \rm{selected}} \times P \times A} \\
   \rho_{r,x}^{\rm{mock}} &=& \frac{n_{r,x}^{\rm{mock}} + f_{r, \rm{mock}} \times P}{N_r^{\rm{mock}} + f_{r, \rm{mock}} \times P \times A} \\

where :math:`A` is the number of characters (e.g., 20 for amino acids), :math:`P > 0` is a pseudocount that is added to each observed count (specified by ``--pseudocount`` option to :ref:`dms2_fracsurvive` / :ref:`dms2_batch_fracsurvive`), and :math:`f_{r, \rm{selected}}` and :math:`f_{r, \rm{mock}}` are defined so that the pseudocount is scaled up for the library (*mock* or *selected*) with higher depth at site :math:`r`:

.. math::

   f_{r, \rm{selected}} = \max\left[1, \left(\sum_x n_{r,x}^{\rm{selected}}\right) / \left(\sum_x n_{r,x}^{\rm{mock}}\right)\right]


.. math::
   f_{r, \rm{mock}} = \max\left[1, \left(\sum_x n_{r,x}^{\rm{mock}}\right) / \left(\sum_x n_{r,x}^{\rm{selected}}\right)\right].

The reason for scaling the pseudocount by library depth is that the *mock* and *selected* libraries may often be sequenced at different depths.
If the same pseudocount is added to both, then the relative frequencies of mutations will be systematically different than one even if the relative counts for the wildtype and mutant amino acid are the same in both two conditions.
Scaling the pseudocounts by the ratio of depths fixes this problem.
If you are using ``--chartype`` of ``codon_to_aa``, the counts for amino acids are aggregated **before** doing the calculations above.

If the mutation has no effect on whether or not a variant survives the selection, then we expect :math:`\rho_{r,x}^{\rm{selected}} = \rho_{r,x}^{\rm{mock}}`.
On the other hand, if the mutation makes the variant more likely to survive the selection, then we expected :math:`\rho_{r,x}^{\rm{selected}} > \rho_{r,x}^{\rm{mock}}`.
If the **only** variants that survive have the mutation, then we expect :math:`\rho_{r,x}^{\rm{selected}} = 1` and the total fraction of the library :math:`\gamma` that survives to just be equal to the unselected frequency of this mutation, :math:`\gamma = \rho_{r,x}^{\rm{mock}}`.
More generally, the fraction of variants with mutation of :math:`r` to :math:`x` that survive the selection is

.. math::

   F_{r,x} = \frac{\gamma \times \rho_{r,x}^{\rm{selected}}}{\rho_{r,x}^{\rm{mock}}}.

It is this quantity :math:`F_{r,x}` that is reported by :ref:`dms2_fracsurvive` / :ref:`dms2_batch_fracsurvive` as the fraction of variants with the mutation surviving the selection.

Error correction
++++++++++++++++++++++++
You can optionally correct for the potential inflation of some counts by sequencing errors by using the ``--err`` option to :ref:`dms2_fracsurvive` / :ref:`dms2_batch_fracsurvive`.
To do this, you should use some control where the **only** mutations will come from sequencing / library-prep-induced errors that are presumed to occur at the same rate as the actual samples being analyzed.
Typically, this would be deep sequencing of a wildtype control.

With an error control, the counts are adjusted as follows.
Let :math:`n_{r,x}` be the counts for character :math:`x` at site :math:`r` in the *mock* or *selected* sample.
Let :math:`n^{\rm{err}}_{r,x}` be the number of counts of :math:`x` at site :math:`r` in the error control.
Define

.. math::

   \epsilon_{r,x} = \left(n^{\rm{err}}_{r,x}\right) / \left(\sum_y n^{\rm{err}}_{r,y}\right).

When :math:`x \ne \operatorname{wt}\left(r\right)` then :math:`\epsilon_{r,x}` is the rate of errors to :math:`x` at site :math:`r`, and when :math:`x = \operatorname{wt}\left(r\right)` then :math:`\epsilon_{r,x}` is one minus the rate of errors away from the wildtype at site :math:`r`.

We then adjust observed counts :math:`n_{r,x}` to the error-corrected counts :math:`\hat{n}_{r,x}` by

.. math::

   \hat{n}_{r,x} = \begin{cases}
   \max\left[\left(\sum_y n_{r,y}\right) \left(\frac{n_{r,x}}{\sum_y n_{r,y}} - \epsilon_{r,x}\right), 0\right] & \mbox{if } x \ne \operatorname{wt}\left(r\right) \\
   n_{r,x} / \epsilon_{r,x} & \mbox{if } x = \operatorname{wt}\left(r\right) \\
   \end{cases}

These adjusted counts are then used in place of the un-adjusted counts in Equation :eq:`fracsurvive_rho` above.

Programs and examples
------------------------------------------
The `dms_tools2`_ software contains programs to calculate the fraction of each mutation surviving:

    * :ref:`dms2_fracsurvive` estimates the fraction surviving from counts files and the overall library fractiont that survives.

    * :ref:`dms2_batch_fracsurvive` runs :ref:`dms2_fracsurvive` on a set of samples and summarizes the results. It is generally easier to use :ref:`dms2_batch_fracsurvive` rather than runnning :ref:`dms2_fracsurvive` directly.

    * :ref:`dms2_logoplot` can be used to create logo plots that visualize the fraction surviving.


.. include:: weblinks.txt

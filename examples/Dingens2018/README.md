# Complete functional mapping of infection- and vaccine-elicited antibodies against the fusion peptide of HIV
This repository analyzes the mutational antigenic profiling data for a study of antibodies against HIV.

The citation for the final published study is [Dingens et al (2018)](https://doi.org/10.1371/journal.ppat.1007159).

The mutational antigenic profiling experiments were performed by Adam Dingens in the [Bloom lab](http://research.fhcrc.org/bloom/en.html) and [Overbaugh lab](https://research.fhcrc.org/overbaugh/en.html) in the summer/autumn of 2017. 

Briefly, we used mutational antigenic profiling to comprehensively map escape from *N123-VRC34.01*, *2712-vFP16.02*, and *2716-vFP20.01* (abbreviated *VRC34*, *FP16-02*, and *FP20-01*) using **BG505.T332N** mutant Env virus libraries.
The generation and characterization of the BG505 mutant virus libraires is described in detail in [Haddox, Dingens et al. eLife 2018](https://elifesciences.org/articles/34420).
The basic mutational antigenic profiling approach for HIV is described in [Dingens et al. Cell Host & Microbe 2017](http://dx.doi.org/10.1016/j.chom.2017.05.003).

The antibodies used here are described in [Xu et al, bioRxiv, DOI 10.1101/306282](https://doi.org/10.1101/306282).

Here we use [dms_tools2](https://jbloomlab.github.io/dms_tools2/) to analyze the data. 
Specifically, the [Jupyter notebook](http://jupyter.org/) [analysis_notebook.ipynb](analysis_notebook.ipynb) does the complete analysis of the mutational antigenic profiling data. 
It downloads the deep sequencing data from the [Sequence Read Archive](http://www.ncbi.nlm.nih.gov/sra), processes this  data, and then analyzes the selection in the context of each antibody. 

The general strategy for barcoded subamplicon Illumina sequencing approach is described [here](https://jbloomlab.github.io/dms_tools2/bcsubamp.html), and the computation of *differential selection* is described [here](https://jbloomlab.github.io/dms_tools2/diffsel.html).

## Organization
The mutational antigenic profiling analysis is performed by the [Jupyter notebook](http://jupyter.org/) [analysis_notebook.ipynb](analysis_notebook.ipynb). 

All the input data needed to run [`analysis_notebook.ipynb`](analysis_notebook.ipynb)) are in the [./data/](./data/). 
These files are described in the iPython notebok when used. 

All results are placed in a created directory named `./results/`.



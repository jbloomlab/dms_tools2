# Deep mutational scanning of the Perth/2009 H3N2 HA

## Overview

This is the analysis code directory for the deep mutational scanning of the A/Perth/16/2009(H3N2) HA. 
We analyze the deep mutational scanning data, and then see if it useful for identifying the fates of viral lineages in nature.

The citation for this study is [Lee et al, 2018](https://doi.org/10.1101/298364).

## Organization

#### Analysis notebooks
There are two analysis notebooks:
  * [analysis_notebook.ipynb](analysis_notebook.ipynb): in this notebook, we analyze deep sequencing data and estimate the amino-acid preferences using [dms_tools2](https://jbloomlab.github.io/dms_tools2/). We also compare the shifts in amino-acid preferences between the Perth/2009 H3 and the WSN/1933 H1 HA's.
  * [max_frequency_analysis.ipynb](max_frequency_analysis.ipynb): in this notebook, we use the Perth/2009 H3 HA and WSN/1933 H1 HA amino-acid preferences to examine the relationships between the frequencies achieved by mutations in a human H3N2 phylogeny and the effects of these mutations.

#### Input files
The [./data/](./data/) subdirectory contains the following input files:

  * `1RVX_trimer_sequentialnumbering.dssp`: the `DSSP` file for the H1 HA (`PDB 1RVX`)
  * `4O5N_trimer.dssp`: the `DSSP` file for the H3 HA (`PDB 4O5N`)
  * `abs_conserved_sites.txt`: list of sites that are absolutely conserved across all 18 HA subtypes
  * `clade_specific_sites.txt`: list of sites that have clade-specific identities in the clade containing H1 and the clade containing H3
  * `domains.csv`: csv file listing each site and its associated HA domain
  * `flu_h3n2_ha_1968_2018_6v_frequencies.json.gz`: `JSON` file of the frequency trajectories of all nodes in a human H3N2 phylogeny (see file below)
  * `flu_h3n2_ha_1968_2018_6v_tree.json.gz`: tree build of human H3N2 influenza virus HA sequences from 1968 to 2018 in a `JSON` file format
  * `flu_h3n2_ha_1968_2018_30v_unpassaged_frequencies.json.gz`: `JSON` file of the frequency trajectories of all nodes in a human H3N2 phylogeny using only unpassaged viruses (see file below)
  * `flu_h3n2_ha_1968_2018_30v_unpassaged_tree.json.gz`: tree build of only unpassaged human H3N2 influenza virus HA sequences from 1968 to 2018 in a `JSON` file format
  * `H1toH3_renumber.csv`: file to convert from H1 sequential numbering to H3 numbering
  * `H3_human_alignment.fa`: file of subsampled human seasonal H3N2 HA sequences for running `phydms`
  * `H3_swine_alignment.fa`: file of subsampled swine H3N2 HA sequences for running `phydms`
  * `H3renumbering_scheme.csv`: file to convert from H3 sequential numbering to H3 numbering
  * `Perth09_HA_reference.fa`: the reference sequence for the wildtype Perth/2009 HA
  * `Perth2009_compareprefs_renumber.csv`: file to renumber the Perth/2009 H3 preferences to allow for comparing the preferences between H3 and H1
  * `Perth2009_WSN_aa_align.fa`: the alignment file for the Perth/2009 H3 and WSN/1933 H1, aligned using `MAFFT`
  * `WSN_avgprefs_rescaled_H3numbering.csv`: the across-replicate and re-scaled WSN/1933 H1 HA amino-acid preferences, in H3 numbering
  * `WSN_compareprefs_renumber.csv`: file to renumber the WSN/1933 H1 preferences to allow for comparing the preferences between H1 and H3
  * `WSN_HA_reference.fa`: the reference sequence for the wildtype WSN/1933 H1 HA
  * `WSN_replicate-*_prefs.csv`: the unscaled amino-acid preferences files for the WSN/1933 H1 HA for each replicate, in sequential numbering
  * `Env_BF520_replicate-*_prefs_rescaled.csv`: amino-acid preferences files for HIV Env (BF520) for each replicate. These amino-acid preferences have already been re-scaled. Taken from [Haddox et al](https://doi.org/10.1101/235630).
  * `BG505_to_BF520_prefs_dist.csv`: shifts between BG505 and BF520 preferences as calculated by [Haddox et al](https://doi.org/10.1101/235630).

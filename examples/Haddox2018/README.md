# Mapping mutational effects along the evolutionary landscape of HIV envelope.

Analysis for "Mapping mutational effects along the evolutionary landscape of HIV envelope", which was published in _eLife_ at [DOI 10.7554/eLife.34420](https://doi.org/10.7554/eLife.34420).

This repository is for the project by Hugh Haddox, Adam Dingens, and [Jesse Bloom](https://research.fhcrc.org/bloom/en.html) analyzing the effects of mutations to HIV Env in different strain backgrounds.
It involved deep mutational scanning of the BG505 and BF520 strains.
Sarah Hilton and Jesse Bloom developed the capabilities to use gamma-distributed models in [phydms](http://jbloomlab.github.io/phydms/).

## Analysis
The actual analysis is performed by the Jupyter notebook [analysis_notebook.ipynb](analysis_notebook.ipynb).
The analysis primarily uses [dms_tools2](https://jbloomlab.github.io/dms_tools2/), and also makes some use of [phydms](http://jbloomlab.github.io/phydms/).

## Input data
All input data required by [analysis_notebook.ipynb](analysis_notebook.ipynb) is in the [./data/](./data/) subdirectory.
Specifically, this subdirectory includes the following files:

  * The wildtype *env* coding sequences for the BG505 and BF520 strains used in the experiments are in [./data/BG505_env.fasta](./data/BG505_env.fasta) and [./data/BF520_env.fasta](./data/BF520_env.fasta).

  * A protein alignment of the Env homologs used in this study plus HXB2 is in [./data/Env_protalignment_manualtweaks.fasta](./data/Env_protalignment_manualtweaks.fasta). This is a manually tweaked version of an alignment created by [mafft](https://mafft.cbrc.jp/alignment/software/). The alignment is used to get the other homologs into the HXB2 numbering scheme. The manual tweaking was done by Hugh Haddox in regions of the variable loops, which are hard to align due to low identity and many indels. Specifically,. Hugh notes that he:

    - re-aligned BF520 in the regions between 184-191 (noninclusive bounds in variable loop 1; HXB2 numbering) and 395-412 (noninclusive bounds in variable loop 4; HXB2 numbering).

    - re-aligned LAI in the region between 137-144 (noninclusive bounds in variable loop 1; HXB2 numbering).


  * An alignment of HIV Env coding sequences is in [./data/HIV1_FLT_2016_env_DNA.fasta](./data/HIV1_FLT_2016_env_DNA.fasta). This alignment is used for the phylogenetic analyses. This alignment was downloaded from the [Los Alamos (LANL) HIV sequence database](http://www.hiv.lanl.gov/). Specifically, it was downloaded using the following settings:

        - Alignment type: Filtered web
        - Organism: HIV-1/SIVcpz
        - Region: Env
        - Subtype: M group without recombinants (A-K)
        - DNA/Protein: DNA
        - Year: 2016
        - Format: FASTA

  * A breakdown of Env into regions (*gp41*, *gp120 variable loops*, *gp120 other regions*) in HXB2 nubmering is in [./data/Env_regions.csv](./data/Env_regions.csv). The variable loop definitions are the [ones provided by the Los Alamos database](https://www.hiv.lanl.gov/content/sequence/VAR_REG_CHAR/variable_region_characterization_explanation.html). The gp41 definition is also the one provided at that site.

  * Files used for calculating the solvent accessibility of Env using [dssp](http://swift.cmbi.ru.nl/gv/dssp/). We calculated solvent accessibility on two Env trimers: [PDB 5FYL](http://www.rcsb.org/pdb/explore.do?structureId=5fyl), which is the "closed" pre-fusion BG505 trimer, and [PDB 5VN3](https://www.rcsb.org/structure/5vn3), which is stabilized in the CD4-bound state. Specifically:

    1. Both structures were downloaded from the PDB.

    2. Next, both PDB files were modified by removing any non-Env chains (i.e., ones from antibodies or CD4). For 5VN3, a complete trimer remained, which we saved in the file [./data/5VN3_Env_trimer.pdb](./data/5VN3_Env_trimer.pdb). For 5FYL, only a single Env monomer remained, which we used to recreate a trimer, as described in the next step.

    3. Next, for 5FYL, we generated symmetry partners using the [PyMol symexp command](https://pymolwiki.org/index.php/Symexp) to create the full Env trimer, re-named the gp120 and gp41 chains so that they are all unique, and removed lines that start with `TER`. The result of all of these operations is in the file [./data/5FYL_Env_trimer_rmTER.pdb](./data/5FYL_Env_trimer_rmTER.pdb).

    3. [dssp](http://swift.cmbi.ru.nl/gv/dssp/) was run on both [./data/5FYL_Env_trimer_rmTER.pdb](./data/5FYL_Env_trimer_rmTER.pdb) and [./data/5VN3_Env_trimer.pdb](./data/5VN3_Env_trimer.pdb) to compute the absolute solvent accessibility of each residue in each structure, and the results were saved to [./data/5FYL_dssp.txt](./data/5FYL_dssp.txt) and [./data/5VN3_dssp.txt](./data/5VN3_dssp.txt), respectively.

  * Files giving the amino-acid preferences for influenza HA (A/WSN/1933 H1N1 strain) to be used for a control comparison to the preferences measured for Env. These are from the experiments of [Doud and Bloom (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27271655), although the files are actually taken from the re-analysis of the data in the [dms_tools2 example](https://github.com/jbloomlab/dms_tools2/blob/master/examples/Doud2016/analysis_notebook.ipynb). The (**not** yet re-scaled by a stringency parameter) preferences are in the following files:

    - [data/Doud2016_HA_replicate-1_prefs.csv](data/Doud2016_HA_replicate-1_prefs.csv)

    - [data/Doud2016_HA_replicate-2_prefs.csv](data/Doud2016_HA_replicate-2_prefs.csv)

    - [data/Doud2016_HA_replicate-3_prefs.csv](data/Doud2016_HA_replicate-3_prefs.csv)
